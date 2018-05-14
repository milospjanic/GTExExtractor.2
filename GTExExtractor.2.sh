#!/bin/bash

FILE="./GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct"
FILE2="./GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
GENES="./genes.txt"
TISSUES="./tissues.txt"

#check if GTEx file exists, if not, download 

if [ ! -f $FILE ]
then
wget https://www.dropbox.com/s/c47ywdnnbge0y4j/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz
gunzip GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz
fi

#check if the Sample Attribution file exists if not
#download GTEx Analysis V6p release, de-identified, open access version of the sample annotations, available in dbGaP.	

if [ ! -f $FILE2 ]
then
wget https://www.dropbox.com/s/23c3zs9igediuy1/GTEx_Data_V6_Annotations_SampleAttributesDS.txt
fi


#grep your tissue of interest, need exact tissue name from the GTEx file
#par=$1 #pass the first argument to a variable so that it can be used in grep with double quotes, in order to use spaced argument 
#grep "$par" GTEx_Data_V6_Annotations_SampleAttributesDS.txt | cut -f1 > $par.sample_IDs

while
IFS= read -r line
do 
echo $line
echo ${line//[[:blank:]]/}
line2=$(echo ${line//[[:blank:]]/})
grep "$line" GTEx_Data_V6_Annotations_SampleAttributesDS.txt | cut -f1 > "$line".sample_IDs
done < "$TISSUES"

ls -1 *.sample_IDs > sample_IDs_all
input="./sample_IDs_all"

#write an extraction script
echo "
awk -v COLT=\$1 '
        NR==1 {
                for (i=1; i<=NF; i++) {
                        if (\$i==COLT) {
                                title=i;
                                print \$i;
                        }
                }
        }
        NR>1 {
                if (i=title) {
                        print \$i;
                }
        }
' \$2
" > extractor.sh
chmod 755 extractor.sh

#create header
cat GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct | head -n3 | tail -n1 > header

#exctract data for each gene in genes.txt from the GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct (all GTEX individuals) 
while
IFS= read -r line
do 
echo $line
grep -P "\t$line\t" GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct > $line.gene.rpkm.txt
cat header $line.gene.rpkm.txt > $line.gene.rpkm.header.txt
done < "$GENES"

#list the file names created in previous step to a file
ls -1 *gene.rpkm.header.txt > gene.file.names.txt
input2="./gene.file.names.txt"

#
while 
IFS= read -r line 
do 
	echo $line >$line.output

		while
		IFS= read -r tissues
		do
        		while 
        		IFS= read -r samples_ID
	
        		do
        		./extractor.sh $samples_ID $line | awk 'NR%2==0'
       			done < "$tissues" >> $line."$tissues".output
			echo "$tissues" | cat - $line."$tissues".output > temp && mv temp $line."$tissues".output
	
		done < "$input"

		paste $line.*.output > $line.results.txt
		sed -i 's/.sample_IDs//g' $line.results.txt

done < "$input2"

#paste *output > results.txt

#sed -i 's/.gene.rpkm.header.txt//g' results.txt

j=1
while
IFS= read -r line
do 
eval "k$j=$line" 
TMP="k$j"
    echo ${!TMP}

j=$(expr $j + 1)


echo "#!/usr/bin/Rscript
library(reshape2)
library (ggplot2)


p<-read.table(\"${!TMP}.gene.rpkm.header.txt.results.txt\", header=T, fill = TRUE, check.names=FALSE, sep=\"\t\")
pmelt<-melt(p)
pmelt2<-na.omit(pmelt)
pdf(\"${!TMP}.output_gtexex.pdf\", width = 7+length(p)/2)
p <- ggplot(pmelt2, aes(x=variable, y=value)) + geom_violin(aes(fill=variable), scale=\"width\") + geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=.5-length(p)/80, binwidth=.7-length(p)/20)+labs(title=\"Expression of ${!TMP} in GTEx samples\", x=\"GTEx tissues\", y=\"RPKM expression level\") + theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1),axis.text.y = element_text(size=16), axis.title.x = element_text(size=18),axis.title.y = element_text(size=18), plot.title=element_text(size=16))
p
dev.off()

" > script.${!TMP}.r

done < "$GENES"

#run R script

find . -name "*script.*.r" | xargs -I % sh -c 'chmod 755 %; %;'

#myFiles <- list.files(pattern="*.results.txt")
#for (k in 1:length(myFiles)){
#}

