bam_dir="/chemin/to/BAM"
output_dir="/chemin/to/save/fastq"

mkdir -p "$output_dir"

for bam in "$bam_dir"/*.bam; do
    sample=$(basename "$bam" .bam)
    echo "Conversion rapide : $sample"

    samtools sort -n -@ 8 "$bam" | \
    samtools fastq -@ 8 \
        -1 "$output_dir/${sample}_R1.fastq.gz" \
        -2 "$output_dir/${sample}_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n
done




for file in /home/pape/Bureau/BAM_bai/*.bam; do
    samtools index "$file"
done

