# GTExExtractor.2

Previously, I have generated a script, GTExExtractor, to show the distribution of expression of multiple genes in a single GTEx tissue in a form of violin plots. Here, I generated a script, GTExExtractor.2, that will show the distribution of expression of a single gene in a multiple GTEx tissue that are selected by the user, and the script will automate this process for a list of input genes, outputing a series of pdfs for each gene.

GTEXExtractor.2 is a combined bash/R script to extract individual level data from the GTEx database, and plot RPKM distributions for the genes of interest in a form of violin plots. GTExExtractor.2 will download individual level data for all GTEx tissues stored in a file **GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct**, in case the file is not present in the working directory. In case the file is present the script will skip the download. The script will read gene names from the genes.txt file provided and extract the RPKM values from all GTEx samples for the provided genes of interest. Values for each gene will be stored in separate files. Tissues of interest have to be stored in a file tissues.txt. Next, sample IDs for the tissues of interest have to be determined. To do that GTExExtractor will check for the file **GTEx_Data_V6_Annotations_SampleAttributesDS.txt**,that contains sample IDs for each GTEx tissue and download it if not already present. The script will extract the sample IDs from the file GTEx_Data_V6_Annotations_SampleAttributesDS.txt for each tissue of interest that was provided in the file tissues.txt. Then, these sample IDs will be used to match the IDs from the individual gene files with RPKM values generated in the previous step, producing tissue-specific expression data for each gene. The files are then combined into a table. The table is the imported into R with Rscript and plotted as a visually representative violin plot.


# Example of usage

Content of the genes.txt

<pre>
SMAD3
AHR
TCF21
</pre>

To run the script 

<pre>
wget https://raw.githubusercontent.com/milospjanic/GTExExtractor/master/GTExExtractor.sh
chmod 755 GTExExtractor.sh
./GTExExtractor.sh 'Artery - Coronary'
</pre>

Check the output pdf file:

![alt text](https://github.com/milospjanic/GTExExtractor/blob/master/output_gtexex.21.png)

If you have more genes in the input file the script will automatically widen the graph

![alt text](https://github.com/milospjanic/GTExExtractor/blob/master/output_gtexex.18.png)

For example, if we run the script for one highly expressed gene like GAPDH and all low expressed TFs the output will look narrower:

![alt text](https://github.com/milospjanic/GTExExtractor/blob/master/output_gtexex.19.png)
