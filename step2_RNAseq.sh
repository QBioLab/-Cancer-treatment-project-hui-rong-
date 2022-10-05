#!/bin/bash

reference=../rawdata/reference/WHR_P185_EG
featureCounts=/datahub/CR/Software/subread-2.0.1-Linux-x86_64/bin/featureCounts
gtf=../rawdata/reference/genes.gtf

pre_mapping()
{
    fastp -w 4 -c -h ${metadata}/${1}.fastp.html -j ${metadata}/${1}.fastp.json \
    -i ${rawdata}/${1}.R1.fq.gz -o ${metadata}/${1}_r1.fq.gz \
    -I ${rawdata}/${1}.R2.fq.gz -O ${metadata}/${1}_r2.fq.gz 
}

mapping()
{
    STAR --runThreadN 10 --genomeDir ${reference} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ${metadata}/${1}_r1.fq.gz ${metadata}/${1}_r2.fq.gz --outFileNamePrefix ${metadata}/${1}.
    samtools index ${metadata}/${1}.Aligned.sortedByCoord.out.bam
    samtools sort -n -o ${metadata}/${1}_nsort.bam ${metadata}/${1}.Aligned.sortedByCoord.out.bam
}

count_gene()
{   
    start_time=$(date +%s)
    ${featureCounts} -a ${gtf} -o ${metadata}/${1}.fc.txt --extraAttributes 'gene_name' -p -T 10 ${metadata}/${1}_nsort.bam
    end_time=$(date +%s)
    cost_time=$[ $end_time-$start_time ]
    echo "featureCounts time is $(($cost_time/60))min $(($cost_time%60))s"
    start_time=$(date +%s)
    htseq-count -n 10 -r pos --additional-attr=gene_name ${metadata}/${1}.Aligned.sortedByCoord.out.bam ${gtf} > ${metadata}/${1}.hc.txt
    end_time=$(date +%s)
    cost_time=$[ $end_time-$start_time ]
    echo "htseq-count time is $(($cost_time/60))min $(($cost_time%60))s"
}

quality_report()
{
    qualimap bamqc --java-mem-size=10G -bam ${metadata}/${1}.Aligned.sortedByCoord.out.bam -ip -nt 10 -outdir ${metadata}/${1}_bamqc 
    #qualimap rnaseq --java-mem-size=10G -bam ${metadata}/${1}_nsort.bam -gtf ${gtf} -pe -s -outdir ${metadata}/${1}_rnaseq
}

rawdata=../rawdata
metadata=../metadata
result=../result

for prefix in Fluc_D1 Fluc_D2 Fluc_D3 Fluc_L1; do 
{
    time pre_mapping "${prefix}" > "${result}/${prefix}.pre_mapping.log" 2>&1
    time mapping "${prefix}" > "${result}/${prefix}.mapping.log" 2>&1
    time count_gene "${prefix}" > "${result}/${prefix}.count_gene.log" 2>&1
    time quality_report "${prefix}" > "${result}/${prefix}.quality_report.log" 2>&1
} & done

for prefix in Fluc_L2 Fluc_L3 mINF_D1 mINF_D3; do 
{
    time pre_mapping "${prefix}" > "${result}/${prefix}.pre_mapping.log" 2>&1
    time mapping "${prefix}" > "${result}/${prefix}.mapping.log" 2>&1
    time count_gene "${prefix}" > "${result}/${prefix}.count_gene.log" 2>&1
    time quality_report "${prefix}" > "${result}/${prefix}.quality_report.log" 2>&1
} & done

for prefix in mINF_D4 mINF_L2 mINF_L3 mINF_L4; do 
{
    time pre_mapping "${prefix}" > "${result}/${prefix}.pre_mapping.log" 2>&1
    time mapping "${prefix}" > "${result}/${prefix}.mapping.log" 2>&1
    time count_gene "${prefix}" > "${result}/${prefix}.count_gene.log" 2>&1
    time quality_report "${prefix}" > "${result}/${prefix}.quality_report.log" 2>&1
} & done