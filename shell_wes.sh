#!/bin/bash

# Read CSV file
file="samplesheet.csv"
tail -n +2 "$file" | while IFS= read -r line; 
do
  # Split line by comma
  IFS=',' read -ra split <<< "$line"
  
  # Store first and second elements into separate variables
  sampleid="${split[0]}"
  read1="${split[1]}"
  read2="${split[2]}"
  
  # remove if the directory already exists
  if [ -d $sampleid ]; then
     rm -rf $sampleid/
  fi
  #make the required directory
  mkdir $sampleid
  cp $sampleid*.fastq.gz "$sampleid"_1.fastq.gz
  cp $sampleid*.fastq.gz "$sampleid"_2.fastq.gz
  cd $sampleid/
  
  # Print the name of sample undergoing analysis 
  echo "Analysis started for $sampleid"

  echo "Quality check $sampleid "
  #fastqc "$read1"
  #fastqc "$read2"
  echo "$sampleid"_1_paired.fastq.gz

  echo "Trimming $sampleid "
  java -jar /Users/ankitverma/Documents/tutorial/dollar_education/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar PE "$read1" "$read2" "$sampleid"_1_paired.fastq.gz  "$sampleid"_1_unpaired.fastq.gz "$sampleid"_2_paired.fastq.gz "$sampleid"_2_unpaired.fastq.gz ILLUMINACLIP:/Users/ankitverma/Documents/tutorial/dollar_education/softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

  echo "Alignment $sampleid "
  /Users/ankitverma/Documents/tutorial/dollar_education/softwares/bwa/bwa mem -t 4 -M ./../genome/ref.fa "$sampleid"_1_paired.fastq.gz "$sampleid"_2_paired.fastq.gz > "$sampleid".sam
  
  echo "Convert SAM to BAM $sampleid"
  /usr/local/bin/samtools view -bS "$sampleid".sam > "$sampleid".bam

  echo "Sort the BAM $sampleid"
  /usr/local/bin/samtools sort -o "$sampleid".sorted.bam "$sampleid".bam

  echo "Index the BAM $sampleid"
  /usr/local/bin/samtools index "$sampleid".sorted.bam

  echo "Remove duplicates $sampleid"
  /opt/homebrew/Cellar/picard-tools/3.0.0/bin/picard MarkDuplicates I="$sampleid".sorted.bam O="$sampleid".sorted.dedup.bam METRICS_FILE="$sampleid"_metrics.txt VALIDATION_STRINGENCY=LENIENT
  /usr/local/bin/samtools index "$sampleid".sorted.dedup.bam


  echo "Variant calling $sampleid"
  /Users/ankitverma/Documents/tutorial/dollar_education/softwares/gatk-4.4.0.0/gatk HaplotypeCaller -R ./../genome/ref.fa -I "$sampleid".sorted.dedup.bam -O "$sampleid".vcf.gz -ERC GVCF

  echo "Variant filtering $sampleid"
  /Users/ankitverma/Documents/tutorial/dollar_education/softwares/gatk-4.4.0.0/gatk SelectVariants -R ./../genome/ref.fa -V "$sampleid".vcf.gz -O "$sampleid".filtered.vcf.gz --select-type-to-include SNP vcftools --gzvcf "$sampleid".filtered.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out "$sampleid".filtered bgzip "$sampleid".filtered.recode.vcf tabix -p vcf "$sampleid".filtered.recode.vcf.gz


  cd ..

done

