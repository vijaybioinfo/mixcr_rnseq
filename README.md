## MiXCR analysis

#### These scripts help you utilise MiXCR on RNA-seq data.

All the information goes into the configuration file (YAML format). There is an example (config.yaml) with comments regarding the files' format.

The files you need to prepare are:
  - You need to provide the path to the location of the FASTQ paths. If not all of them need to be analysed, you just need to prepare a text (CSV, TXT) file with paths to the FASTQ files.

### Install
Clone this repository (your ~/bin folder is a good place).
```
git clone https://github.com/vijaybioinfo/mixcr_rnseq.git
cd mixcr_rnseq
```

Make sure your config template is pointing to where you have the pipeline.
This will also add the run.sh script as an alias.
```
sh locate_pipeline.sh
```

### Run the pipeline
After you've added the necessary information to the YAML file you can call the pipeline.
```
mixcr_rnseq -y /path/to/project/config.yaml
