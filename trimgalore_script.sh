mkdir fastqc_results

for file in ./fastq_files/*
do 
	trim_galore \
		--quality 25 \
		--stringency 5 \
		--length 50 \
		--output_dir ./fastqc_results \
		--fastqc \    # output dir will contains both the trimming logs and fastqc results
		${file}
done
