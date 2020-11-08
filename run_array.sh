#!/usr/bin/R

#########
# MIXCR #
#########

# set -euo pipefail

function usage () {
    cat >&2 <<EOF

USAGE: $0 [-y] [options]
  -y <config file> : Path to the YAML config file. Required.
  -v Verbose.
  -h Print the usage info.

EOF
}
# initial : makes this loop silent and now requires '?)'
# ${opt} is each option and ${OPTARG} its the argumet (if a colon is there ${opt}:)
VERBOSE=FALSE
while getopts ":y:vh" opt; do
  case ${opt} in
    y) CONFIG_FILE=${OPTARG};;
    v) VERBOSE=TRUE;;
    h) usage; exit 1;;
    \?) echo "No -${OPTARG} argument found."; usage; exit 1;;
  esac
done
if [[ ${OPTIND} -eq 1 ]] ; then
    usage; exit 1
fi

#### Parameters #### -----------------------------------------------------------
function read_yaml(){
  sed 's/#.*//g' ${1} | grep ${2}: | sed 's/.*:[^:\/\/]//; s/\"//g'
}
FASTQS="$(read_yaml ${CONFIG_FILE} fastqs)"
OUTPUT_DIR="$(read_yaml ${CONFIG_FILE} output_dir)"
PROJECT_ID="$(read_yaml ${CONFIG_FILE} project_id)"
OUTDIR=${OUTDIR%/}/${PROJECT_ID}

## Job coordination
SUBMIT="$(read_yaml ${CONFIG_FILE} submit)"
DEPEND="$(read_yaml ${CONFIG_FILE} dependency)"
WALLTIME="$(read_yaml ${CONFIG_FILE} walltime)"
MEM="$(read_yaml ${CONFIG_FILE} mem)"
NODES="$(read_yaml ${CONFIG_FILE} nodes)"
PPN="$(read_yaml ${CONFIG_FILE} ppn)"

## Digeting parameters
DEPEND="${DEPEND:-NONE}"

echo; echo "**** Vijay Lab - LJI"
echo -e "\033[0;36m------------------------------- PRESENTING PARAMETERS -------------------------------\033[0m"
echo "Fastqs: ${FASTQS}"
echo "Output directory: ${OUTDIR}"
echo -e "\033[0;36m------------------------------- --------------------- -------------------------------\033[0m"
if [ ${DEPEND} == "NONE"  ]; then read -n 1 -s -r -p "Press any key to continue"; fi; echo
echo; echo
if [[ ! -d "${OUTDIR}" ]]; then mkdir --parents "${OUTDIR}"; fi
ls -loh "${OUTDIR}"
echo -e "\033[0;36m------------------------------- Initialising MIXCR analysis -------------------------------\033[0m"

# Extracting sample names
if [[ -d "${FASTQS}" ]]; then
  echo "From directory"
  FASTQS=${FASTQS%/}
  MYSAMPLES=(`find "${FASTQS}" | grep 001.fastq | grep -v Undetermined | sed 's/_S[0-9]\{1,\}_L.*//g' | sort -u`)
elif [[ ""${FASTQS}"" =~ ".csv" ]]; then
  echo "From file"
  MYSAMPLES=(`cut -d, -f1 "${FASTQS}"`)
else
  MYSAMPLES=${FASTQS[@]}
fi
SAMPLES_FILE="${OUTDIR}/sample_names_${#MYSAMPLES[@]}.txt"
if [[ -s "${SAMPLES_FILE}" ]]; then rm "${SAMPLES_FILE}"; fi # delete it if it exists
echo -e "\033[0;35m>>>>>>>\033[0m Getting samples \033[0;35m<<<<<<<\033[0m"
for SAMPLE_I in ${MYSAMPLES[@]}; do # check which samples need to be processed
  SAMPLE_NAME=$(basename ${SAMPLE_I})
  # echo -e "\033[0;35m>>>>>>>\033[0m ${SAMPLE_NAME} \033[0;35m<<<<<<<\033[0m"
  if [ `ls "${SAMPLE_I}*fastq.gz" | wc -l` -eq 0 ]; then
    echo "${SAMPLE_NAME} does not exist"; continue
  elif [ -s "${OUTDIR}/${SAMPLE_NAME}/clones_TRB.txt" ]; then
    printf "\033[0;32m${SAMPLE_NAME} already processed\033[0m: "; continue
  fi
  echo "${SAMPLE_I}" >> "${SAMPLES_FILE}"
done
SAMPLES2RUN=$(wc -l "${SAMPLES_FILE}" | sed 's/ .*//g')
echo -e "\033[0;35m>>>>>>>\033[0m ${SAMPLES2RUN} ready"

JOBFILE="${OUTDIR}/script_${SAMPLES2RUN}"
cat <<EOT > ${JOBFILE}.sh
#PBS -t 1-${SAMPLES2RUN}%200
#PBS -N MIXCR_${SAMPLES2RUN}
#PBS -o ${JOBFILE}.out.txt
#PBS -e ${JOBFILE}.err.txt
#PBS -m ae
#PBS -M ${USER}@lji.org
#PBS -q default
#PBS -l nodes=${NODES}:ppn=${PPN}
#PBS -l mem=${MEM}
#PBS -l walltime=${WALLTIME}

# This job structure is based of:
# https://learn.lji.org/display/BIODOCS/Creating+a+Job
# https://learn.lji.org/display/BIODOCS/Job+Arrays

######################################################################
#                                                                    #
#   Preface of operations: introduce the code you need to prepare    #
#   your job's parameters; this is useful especially when            #
#   when you have an array                                           #
#                                                                    #
######################################################################

source /home/ciro/scripts/functions/bash_functions.sh
SAMP=\$(head -n \$PBS_ARRAYID "${SAMPLES_FILE}" | tail -n 1)
SNAME=\`basename \${SAMP}\`

echo -----------------------------------------------------------------
echo -n 'Job is running on node '; cat \$PBS_NODEFILE
echo -----------------------------------------------------------------
echo PBS: qsub is running on \$PBS_O_HOST
echo PBS: originating queue is \$PBS_O_QUEUE
echo PBS: executing queue is \$PBS_QUEUE
echo PBS: working directory is \$PBS_O_WORKDIR
echo PBS: execution mode is \$PBS_ENVIRONMENT
echo PBS: job identifier is \$PBS_JOBID
echo PBS: job name is \$PBS_JOBNAME
echo PBS: node file is \$PBS_NODEFILE
echo PBS: current home directory is \$PBS_O_HOME
echo PBS: PATH = \$PBS_O_PATH
echo -----------------------------------------------------------------

######################################################################
#                                                                    #
#   To minimize communications traffic, it is best for your job      #
#   to work with files on the local disk of the compute node.        #
#   Hence, one needs to transfer files from your permanent home      #
#   directory tree to the directory ${WORKDIR} automatically         #
#   created by PBS on the local disk before program execution,       #
#   and to transfer any important output files from the local        #
#   disk back to the permanent home directory tree after program     #
#   execution is completed.                                          #
#                                                                    #
######################################################################

# The working directory for the job is inside the scratch directory
WORKDIR=/mnt/BioScratch/\${USER}/mixcr/\${SNAME}

# This is the directory on 'local' where your project is stored
PROJDIR=${OUTDIR}/\${SNAME}

######################################################################
#                                                                    #
#   Extra job monitoring code.                                       #
#                                                                    #
######################################################################

# Job resource monitoring and summarisation
# https://confluence.lji.org/display/BIODOCS/Job+Resource+Monitoring+and+Summarization
declare -ix RSAMPLER_FREQ_SEC=5
declare -x  RSAMPLER_TMP_PATH=/mnt/BioScratch/\${USER}/mem_usage
declare -x  RSAMPLER_FILE="\${PBS_JOBNAME}_\${PBS_JOBID}"_resource_samples.txt
declare -x  RSAMPLER_OUTPUT="\${RSAMPLER_TMP_PATH}/\${RSAMPLER_FILE}"
declare -x  RSAMPLER_TERM="\${RSAMPLER_TMP_PATH}/\${PBS_JOBID}_finished"
mkdir --parents \${RSAMPLER_TMP_PATH}

# start the sampler starting with our scripts parent id then capture the sampler pid
resourceSampler \$PPID &
RSAMPLER_PID=\$!
RSTART=\$(date +%s)

function logger(){
    echo "$(hostname): Array Number: \$PBS_ARRAYID $(date): \$1"
}

######################################################################
#                                                                    #
#   Transfer files from server to local disk.                        #
#                                                                    #
######################################################################

stagein()
{
   echo ' '
   echo Transferring files from server to compute node
   echo Creating the working directory: \${WORKDIR}
   mkdir --parents \${WORKDIR}
   mkdir --parents \${PROJDIR}
   cd \${WORKDIR}
   echo Writing files in node directory  \`pwd\`

   # Add code here to transfer files from your home directory
   # on the head node to the working directory, e.g:
   #
   # cp -R \${PROJDIR}/directory_to_copy_from ./directory_to_copy_to

   echo Files in node work directory are as follows:
   ls -loh
}

######################################################################
#                                                                    #
#   Execute the run.  Do not run in the background.                  #
#                                                                    #
######################################################################

runprogram()
{
  echo "**** Vijay Lab - LJI"
  echo "Working at: \`pwd\`"
  mixcr='java -Xmx4g -Xms3g -jar /mnt/BioHome/ciro/bin/mixcr-2.1.10/mixcr.jar'

  FILES=(\`ls \${SAMP}*R*.fastq*\`)
  echo "\${#FILES[@]} samples"
  if [ \${#FILES[@]} -gt 2 ]; then
    READSF=(R1 R2)
    for READS_I in \${READSF[@]}; do
      FILES_FOUND=(\`grep \${SAMP} "${SAMPLES_FILE}" | grep \${READS_I}_001.fastq \`)
      cat \${FILES_FOUND[@]} > \${SNAME}.\${READS_I}.merged.fastq.gz
    done
    FFILE1=\${SNAME}.R1.merged.fastq.gz
    FFILE2=\${SNAME}.R2.merged.fastq.gz
  else
    FFILE1=\${FILES[0]}
    FFILE2=\${FILES[1]}
  fi
  echo "Input files:"
  echo "\${FFILE1}"
  echo "\${FFILE2}"

  # The	following lines	are constructed	from the Analysis of RNA-Seq data part of the Quick start
  # https://mixcr.readthedocs.io/en/master/quickstart.html

  echo ">>>>>>>>>>>> Starting MIXCR analysis <<<<<<<<<<<<"
  echo "------------ Aligning"
  \$mixcr align -s hsa -p rna-seq -OallowPartialAlignments=true --verbose \${FFILE1} \${FFILE2} alignments.vdjca
  echo; echo "------------ Assemble"
  \$mixcr assemblePartial alignments.vdjca alignmentsRescued_1.vdjca
  echo; echo "------------ Rescue"
  \$mixcr assemblePartial alignmentsRescued_1.vdjca alignmentsRescued_2.vdjca
  echo; echo "------------ Extended aligments"
  \$mixcr extendAlignments alignmentsRescued_2.vdjca alignmentsRescued_2_extended.vdjca
  echo; echo "------------ Clones"
  \$mixcr assemble alignmentsRescued_2_extended.vdjca clones.clns
  echo; echo "------------ Export"
  \$mixcr exportClones -o --preset full -vHit -jHit -dHit -t clones.clns clones.txt
  echo; echo "------------  - Alpha"
  \$mixcr exportClones -c TRA -o --preset full -vHit -jHit -dHit -t clones.clns clones_TRA.txt
  echo; echo "------------  - Beta"
  \$mixcr exportClones -c TRB -o --preset full -vHit -jHit -dHit -t clones.clns clones_TRB.txt
  echo; echo "DONE"
}

######################################################################
#                                                                    #
#   Copy necessary files back to permanent directory.                #
#                                                                    #
######################################################################

stageout()
{
  echo ' '
  echo Transferring files from compute nodes to server
  echo Writing files in permanent directory \${PROJDIR}

  # if you've created an output directory inside your working
  # directory, you would copy it back to your home directory
  # like so:
  cp -R \${WORKDIR}/clones* \${WORKDIR}/alignments* \${PROJDIR}/

  echo Final files in permanent data directory:
  cd \${PROJDIR}
  ls -loh

  # if [[ condition ]]; then
  #   echo Removing the temporary directory from the compute node
  #   rm -rf \${WORKDIR}
  # fi
}

######################################################################
#                                                                    #
#   The "qdel" command is used to kill a running job.  It first      #
#   sends a SIGTERM signal, then after a delay (specified by the     #
#   "kill_delay" queue attribute (set to 60 seconds), unless         #
#   overridden by the -W option of "qdel"), it sends a SIGKILL       #
#   signal which eradicates the job.  During the time between the    #
#   SIGTERM and SIGKILL signals, the "cleanup" function below is     #
#   run. You should include in this function commands to copy files  #
#   from the local disk back to your home directory.  Note: if you   #
#   need to transfer very large files which make take longer than    #
#   60 seconds, be sure to use the -W option of qdel.                #
#                                                                    #
######################################################################

early()
{
  echo ' '
  echo ' ############ WARNING:  EARLY TERMINATION #############'
  echo ' '
}
trap 'early; stageout' 2 9 15

######################################################################
#                                                                    #
#   Staging in, running the job, and staging out                     #
#   were specified above as functions.  Now                          #
#   call these functions to perform the actual                       #
#   file transfers and program execution.                            #
#                                                                    #
######################################################################

logger "Starting.."
stagein
runprogram
stageout
logger "Finished."

######################################################################
#                                                                    #
#   The epilogue script automatically deletes the directory          #
#   created on the local disk (including all files contained         #
#   therein.                                                         #
#                                                                    #
######################################################################

# stop sampler which will print its info and clean up
touch "\$RSAMPLER_TERM"
# wait for sampler to finish
wait \$RSAMPLER_PID
# rm -f "\$RSAMPLER_OUTPUT"
REND=\$(date +%s)
echo "Start \$(date -d @\${RSTART})"
echo "End \$(date -d @\${REND})"
echo "Finished \${PBS_JOBNAME}: \$PBS_JOBID in \$((REND-RSTART)) seconds."
echo "\`convertsecs \$((REND-RSTART))\` total"

exit
EOT

if [[ -s ${JOBFILE}.out.txt ]]; then rm ${JOBFILE}.*.txt; fi
echo -e "\033[0;36mCheck it out:\033[0m ${JOBFILE}.sh"
if [[ "${SUBMIT}" != "TRUE" ]]; then
  echo "Submit it manually"
else
  if [[ ${DEPEND} != "NONE" ]]; then
    CID=$(qsub -W depend=afterok:${DEPEND} ${JOBFILE}.sh)
  else
    CID=$(qsub ${JOBFILE}.sh)
  fi; CID=`echo ${CID} | sed 's/\..*//'`
  echo ">>>>>> Job ID: ${CID}"
fi
echo -e "\033[0;36m------------------------------- --------------------------- -------------------------------\033[0m"; echo
