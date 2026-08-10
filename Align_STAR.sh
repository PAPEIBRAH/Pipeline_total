#!/bin/bash


########################################
# PARAMETRES
########################################

genome_dir="/home/pape/NAS_1/Pape/Data/Index_gtf50"

fastq_dir="/home/pape/NAS_1/Pape/Data/Fastq_patient_MMG"

out_dir="/home/pape/Bureau/BAM_PATIENTGTF50"


threads=16
jobs=2


mkdir -p "$out_dir"



########################################
# ALIGNEMENT
########################################


align_sample () {


R1="$1"

sample=$(basename "$R1" _R1.fastq.gz)

R2="${fastq_dir}/${sample}_R2.fastq.gz"


if [ ! -f "$R2" ]
then
    echo "FASTQ R2 manquant pour $sample"
    exit 1
fi


echo "================================="
echo "Début : $sample"
date
echo "================================="

# Alignement de STAR optimal pour etudier apres l'epissage 

STAR \
--runThreadN "$threads" \
--genomeDir "$genome_dir" \
--readFilesIn "$R1" "$R2" \
--readFilesCommand zcat \
--twopassMode Basic \
--outBAMcompression 8 \
--outFileNamePrefix "${out_dir}/${sample}_" \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS NM MD \
--outSAMstrandField intronMotif \
--outSJfilterOverhangMin 8 8 8 8 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMultimapNmax 20 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--outSAMunmapped None \
--limitBAMsortRAM 80000000000



# nettoyage fichiers inutiles

rm -f "${out_dir}/${sample}_Log.out"
rm -f "${out_dir}/${sample}_Log.progress.out"



echo "Terminé : $sample"
date


}



export -f align_sample

export genome_dir
export fastq_dir
export out_dir
export threads



########################################
# EXECUTION PARALLELE
########################################


find "$fastq_dir" -name "*_R1.fastq.gz" | \
parallel -j "$jobs" align_sample {}



echo "Tous les alignements sont terminés."
