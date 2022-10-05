#!/bin/bash

cd ../rawdata/reference/
cat vector.fa >> genome.fa
cat vector.gtf >> genes.gtf
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir WHR_P185_EG --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbOverhang 149 