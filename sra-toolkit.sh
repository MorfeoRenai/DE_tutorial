VAR=$(cut -d ',' -f 1 SraRunTable.txt | tail -n +2)

for i in ${VAR}
    do
        echo "Download SRA sample: ${i}"
        prefetch ${i}
        fastq-dump --gzip --defline-qual '+' ${i}
    done
