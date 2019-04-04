#!/bin/bash

# For the STAR alignments we will use the UMItools timmed files in the HiSTA2 directory.
DATADIR='/projects/polysomeProf_BetaCells_Glu/data/reads/'
RESULTDIR='/projects/polysomeProf_BetaCells_Glu/results/201807_STAR/'

## STAR alignments
#for infile in ${DATADIR}*.fastq;
#do
#  FNAME=`echo ${infile} | cut -d'.' -f1`;
#  FNAME=$(basename $FNAME)"/"
#  DIRALIGN=${RESULTDIR}${FNAME}
#  mkdir -p ${DIRALIGN}
#  STAR --runThreadN 11 --genomeDir /projects/STAR_indices/hsapiens/ensembl_38_92_star260 --genomeLoad LoadAndKeep --readFilesIn ${infile} --outFileNamePrefix ${DIRALIGN} --outFilterMultimapNmax 20 --outSAMprimaryFlag AllBestScore --limitBAMsortRAM 32000000000 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outStd Log --seedSearchStartLmaxOverLread 0.5 --winAnchorMultimapNmax 36 --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5 >> STARpipeline.log
#  echo "== Finish run for ${FNAME} ==" >> STARpipeline.log
#done;
# removed parameter --sjdbGTFfile /projects/annotations/GTF/hsapiens/Ensembl_38_92/Homo_sapiens.GRCh38.92.gtf 

exps="heavy_low01 light_high02 light_low03 mono_low01 total_high02 total_low03 heavy_high01 heavy_low02 light_high03 mono_high01 mono_low02 total_high03 heavy_high02 heavy_low03 light_low01 mono_high02 mono_low03 total_low01 heavy_high03 light_high01 light_low02 mono_high03 total_high01 total_low02"

## RSEM analysis.
#cd RSEM_quant
#for exp in ${exps};
#do
#  echo "rsem-calculate-expression --alignments --no-bam-output -p 11 ../${exp}/Aligned.toTranscriptome.out.bam ref/hsapiens_ENSEMBL92 ${exp} >> rsemlog 2>&1"
#  rsem-calculate-expression --alignments --no-bam-output -p 12 ../${exp}/Aligned.toTranscriptome.out.bam ref/hsapiens_ENSEMBL92 ${exp} >> rsemlog 2>&1
#done;

#stringtie to quantify transcripts isoforms
#cd STRINGTIE_quant
#for exp in ${exps};
#do
#  mkdir ${exp};
#  cd ${exp};
#  #echo "stringtie ../../${exp}/Aligned.sortedByCoord.out.bam -o out.gtf -G /projects/annotations/GTF/hsapiens/Ensembl_38_92/Homo_sapiens.GRCh38.92.gtf -e -B -A gene_abund.tab";
#  stringtie ../../${exp}/Aligned.sortedByCoord.out.bam -o out.gtf -G /projects/annotations/GTF/hsapiens/Ensembl_38_92/Homo_sapiens.GRCh38.92.gtf -e -B -A gene_abund.tab -p 12;
#  cd ../;
# done

#remove line containig exon info
cd STRINGTIE_quant
for exp in ${exps};
do
	cd ${exp};
#	echo `pwd`
# 	echo "grep -v "exon_number" out.gtf > transcript.out.gtf"
#	grep -v "exon_number" out.gtf > transcript.out.gtf
cut -f 9 transcript.out.gtf |cut -f 1,2,3,4,6 -d ";" | cut -f 2,4,6,8,10 -d '"' | sed 's/"/\t/g' > transcript_out.tab
	cd ../;
done


# merge counts.
#for infile in ${DATADIR}*genes.results;
#do
#    echo ${infile} | cut -d'.' -f1 >> TPMrowname.tsv;
#done;

# create table for TPM (files are taken in alphabetic order)
#paste *.genes.results | awk 'BEGIN { OFS = " " }{ print $6,$13,$20,$27,$34,$41,$48,$55,$62,$69,$76,$83,$90,$97,$104,$111,$118,$125,$132,$139,$146,$153,$160,$167}' > TPM.txt

## All these crap was from the RIP analysis.

# Alignement with Hisat2
# for infile in ${RESULTDIR}*_trimmed.fastq;
# do
#   FNAME=`echo ${infile} | cut -d'.' -f1`;
#   hisat2 -q -p 12 -5 3 -x /projects/STAR_indices/grch38_Ensembl_HISAT_index/genome -U ${FNAME}.fastq -S ${FNAME}.sam > ${FNAME}.log;
# done;

# SAM BAM conversions, sorting and indexing.
#for infile in ${RESULTDIR}Aligned.sortedByCoord.out.bam;
#do
#  FNAME=`echo ${infile};
#  samtools view -S -b -@ 12 ${FNAME}.sam > ${FNAME}.bam;
#  samtools sort -@ 12 ${FNAME}.bam ${FNAME}_sorted;
#  samtools index ${FNAME}_sorted.bam ${FNAME}_sorted.bai
# done;

## UMITOOLS de-duplication.
#for infile in ${RESULTDIR}_trimmedAligned.sortedByCoord.out.bam;
#do
#   FNAME=`echo ${infile} | cut -d'.' -f1`;
#   echo ${FNAME}
##   samtools index ${FNAME}_sorted.bam ${FNAME}_sorted.bai
##   umi_tools dedup -I ${infile} --output-stats=${FNAME}_stats -S ${FNAME}_deduplicated.bam
#done;

# SAM BAM conversions, sorting and indexing.
# for infile in ${RESULTDIR}*_trimmed.sam;
# do
#   FNAME=`echo ${infile} | cut -d'.' -f1`;
#   samtools view -S -b -@ 12 ${FNAME}.sam > ${FNAME}.bam;
#   samtools sort -@ 12 ${FNAME}.bam ${FNAME}_sorted;
#   samtools index ${FNAME}_sorted.bam ${FNAME}_sorted.bai
# done;


# # BAM sorting and indexing.
# for infile in ${RESULTDIR}*_deduplicated.bam;
# do
#   FNAME=`echo ${infile} | cut -d'.' -f1`;
#   samtools sort -@ 12 ${FNAME}.bam ${FNAME}_sorted;
#   samtools index ${FNAME}_sorted.bam ${FNAME}_sorted.bai;
# done;


