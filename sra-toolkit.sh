VAR=$(cut -d ',' -f 1 SraRunTable.txt | tail -n +2) # select the first field, containg the IDs needed for the sra-toolkit, and eliminate the header "Run"

for i in ${VAR}
    do
        echo "Download SRA sample: ${i}"
        prefetch ${i}                               # for each run ID we prefetch the data...
        fastq-dump --gzip --defline-qual '+' ${i}   # ...and then download it already decrompressed 
    done
    
