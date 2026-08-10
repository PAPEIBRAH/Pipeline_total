# Generation du bed à partir du gtf qui a éte utilise pour l'alignement pour faire Rseqc
conda activate ucsc_tools
which gtfToGenePred
which genePredToBed
cd /home/pape/Bureau
gtfToGenePred \
gencode.v50.annotation.gtf \
gencode.v50.annotation.genePred \
-allErrors
genePredToBed \
gencode.v50.annotation.genePred \
gencode.v50.annotation.bed
head -2 gencode.v50.annotation.bed

#!/bin/bash

BAM_DIR="/home/pape/Bureau/BAM_PATIENTGTF50/BAM"

REF_BED="/home/pape/Bureau/gencode.v50.annotation.bed"

OUTDIR="/home/pape/Bureau/RSeQC_JUNCTION_BAM_PATIENTGTF50"

mkdir -p "$OUTDIR"


for BAM in ${BAM_DIR}/*Aligned.sortedByCoord.out.bam
do

    SAMPLE=$(basename "$BAM" _Aligned.sortedByCoord.out.bam)

    echo "Analyse jonctions : $SAMPLE"


    junction_annotation.py \
    -i "$BAM" \
    -r "$REF_BED" \
    -o "${OUTDIR}/${SAMPLE}"


done


echo "Terminé"








