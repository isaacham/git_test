#Set the working directory
setwd('/Users/admin/Desktop/RNA-seq/tutorial')

#Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)

#Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)



#Data Wrangling

#DESeq expects the counts to have gene IDS as row names
head(counts)
rownames(counts) = counts$Gene.ID
head(counts)

#Remove unused columns (gene ID and gene name)
genes = counts[, c("Gene.ID", "Gene.Name")] #Keeping a copy
counts = counts[, -c(1,2)]
head(counts)

#DESeq epects the metdata matrix to have sample ID's in the rownames
head(metadata)
rownames(metadata) = metadata$Run
head(metadata)

#Only keep columns of interest 
metadata = metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]
metadata

#Rename Sample.Characteristic.genotype to Genotype
colnames(metadata) = c('genotype')
metadata

#Remove spaces in names to avoid DESeq Errors
metadata$genotype[metadata$genotype == 'wild type genotype'] = 'wildtype'
metadata$genotype[metadata$genotype == 'Snai1 knockout'] = 'knockout'
metadata
 
#Turn Genotype into a factor (????) 
metadata$genotype = factor(metadata$genotype, levels=c('wildtype','knockout'))
metadata$genotype



#Spot Checking Expression for knockout gene Snai1
gene_id = genes$Gene.ID[ genes$Gene.Name == 'SNAI1']
gene_counts = counts[gene_id, ]
gene_counts


#Combining metdata columns into the dataframe
gene_data = cbind(metadata, counts=as.numeric(gene_counts))
gene_data


library(ggplot2)
ggplot(gene_data, aes(x = genotype, y=counts, fill=genotype)) + geom_boxplot()

#finished Data Wrangling!!! 



#RUN DESEQ


dds <- DESeqDataSetFromMatrix(countData=counts, colData = metadata, design= ~genotype)

#Ignore gene with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

#Run DESeq
dds <- DESeq(dds)

#Compare expression
res = results(dds, contrast=c('genotype', 'knockout', 'wildtype'), alpha=1e-5)
res





