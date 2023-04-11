#Build databases
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/

annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/ 

annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/

table_annovar.pl example/ex1.avinput humandb/ -buildver hg38 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish  -vcfinput

