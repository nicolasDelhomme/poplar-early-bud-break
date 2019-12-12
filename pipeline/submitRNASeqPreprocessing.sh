#!/bin/bash -l

set -ex # 'set -e' arrête le script quand il y a une erreur à une étape pour éviter de continuer les étapes suivantes pour rien avec des erreurs
		# 'set -x' print et interprete les variables pour se rendre compte là où on a mal tapé un nom de variable

# environment variables
proj=u2017009

mail="katja.stojkovic@slu.se"

in=/mnt/picea/projects/aspseq/u2017009/raw

out=/mnt/picea/projects/aspseq/u2017009

genome=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/STAR/2.5.2b/Potra01

gff3=/mnt/picea/storage/reference/Populus-tremula/v1.1/gff3/Potra01-gene-synthetic-transcripts-wo-intron.gff3

gtf=/mnt/picea/storage/reference/Populus-tremula/v1.1/gtf/Potra01-gene-mRNA-wo-intron.gtf

kallistoFasta=/mnt/picea/storage/reference/Populus-tremula/v1.1/fasta/Potra01-mRNA.fa
			  
kallistoIndex=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/kallisto/Potra01-mRNA.fa.inx

start=9
end=9
mem=256 # minimum 160 Go car c'est la taille du génome indexé pour kallisto qui doit rentrer entièrement en mémoire, donc on prévoit un peu plus de buffer


module load bioinfo-tools samtools fastQvalidator FastQC sortmerna Trimmomatic star htseq kallisto 
# rajouter modules: fastQValidator, FastQC, SortMeRNA, Trimmomatic, STAR, HTSeq-count, Kallisto
#module avail

if [ -z $UPSCb ]; then # -z teste si on a la variable UPSCb dans notre environnement
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir." 
    exit 1
fi

# loop to get every PE fastq file and apply the script runRNASeqPreprocessing.sh to it
for f in `find $in -name "*_[1,2].fastq.gz"`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line; 
# uniq - report or omit repeated lines (With no options, matching lines are merged to the first occurrence.)
do 
bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -m $mem -g $genome -G $gtf -H $gff3 -f $kallistoFasta -K $kallistoIndex -t $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out 
# options of runRNASeqPreprocessing.sh, run it for every PE fastq file called f in the loop
done
