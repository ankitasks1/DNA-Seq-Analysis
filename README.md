# DNA-Seq-Analysis
DNA sequencing analysis

#### The shell script for samplesheet.csv is present in shell_wes.sh
# saved as ./steps_shell_wes.sh

#------------WES/WGS—---------#

# Step 0: Getting data from Public
Get the data from public repository (SRA)
#Install SRA tool kit

#https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

#Select Ubuntu Linux 64 bit architecture/ mac OSX

#Install sra-toolkit

/home/ankits/sratoolkit/bin/fastq-dump --split-files --gzip SRR9876543

# Step 1: Quality Control
fastqc -o fastqc_output/ input.fastq.gz

# Step 2: Trimming
java -jar trimmomatic.jar PE -phred33 input_1.fastq.gz input_2.fastq.gz output_1_paired.fastq.gz output_1_unpaired.fastq.gz output_2_paired.fastq.gz output_2_unpaired.fastq.gz ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3: Mapping
bwa index reference.fa

bwa mem -t 4 -M reference.fa output_1_paired.fastq.gz output_2_paired.fastq.gz > output.sam

# Step 4: Processing SAM/BAM files
samtools view -bS output.sam > output.bam

samtools sort -o output.sorted.bam output.bam

samtools index output.sorted.bam

# Step 5: Removing duplicates

java -jar picard.jar MarkDuplicates I=output.sorted.bam O=output.sorted.dedup.bam METRICS_FILE=metrics.txt VALIDATION_STRINGENCY=LENIENT

samtools index output.sorted.dedup.bam

# Step 6: Variant Calling
#Germline

gatk HaplotypeCaller -R reference.fa -I output.sorted.dedup.bam -O output.vcf.gz -ERC GVCF

#Somatic

gatk Mutect2 -R reference.fasta -I tumor.bam -tumor tumor_sample_name -I normal.bam -normal normal_sample_name --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.vcf.gz -O output.vcf.gz


# Step 7: Variant Filtering
gatk SelectVariants -R reference.fa -V output.vcf.gz -O output.filtered.vcf.gz --select-type-to-include SNP
vcftools --gzvcf output.filtered.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --out output.filtered
bgzip output.filtered.recode.vcf
tabix -p vcf output.filtered.recode.vcf.gz

# Step 8: Annotation
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/

annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/ 

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/

table_annovar.pl sample.avinput humandb/ -buildver hg38 -out myanno_sample -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish  -vcfinput
