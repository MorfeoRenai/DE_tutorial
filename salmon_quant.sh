GEO=$(cat ../geo_accessions.txt)

mkdir ../salmon_results

for i in ${GEO}
do
	SRR=$(grep ${i} ../SraRunTable.txt | cut -d "," -f 1)
	SRR=$(echo ${SRR} | sed "s/ /_trimmed.fq.gz /g")
	SRR=${SRR}_trimmed.fq.gz
	salmon quant -i ../hg38_index --libType A -o ../salmon_results/${i} -r ${SRR}

done
