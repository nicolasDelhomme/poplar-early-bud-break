#!/bin/bash

set -ex

proj=u2017009
mail="katja.stojkovic@slu.se"
in=/mnt/picea/projects/aspseq/$proj/fastqc/raw
out=/mnt/picea/projects/aspseq/$proj/fastqc/raw/multiqc

if [ ! -d $out ]; then
	mkdir -p $out
fi

module load bioinfo-tools multiqc

sbatch --mail-user=$mail -o $in/multiqc.out -e $in/multiqc.err -A $proj $UPSCb/pipeline/runMultiQC.sh $in $out
