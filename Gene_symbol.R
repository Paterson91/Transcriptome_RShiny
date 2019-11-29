#One way of getting EntrezID and Gene Symbols
#PROBLEM: gets rid of lots of genes.
#THIS IS WORKING.Although still can't get it into the heatmap because there are duplicates of the gene symbols?
#Creating accession codes
library(data.table)
library(AnnotationDbi)
library(org.Rn.eg.db)
library("annotate")
library(rae230a.db)

#Set working directory
setwd("/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/GSE8491/")
# load series and platform data from GEO
library(GEOquery)
getGEOSuppFiles("GSE8491")

untar("GSE8491/GSE8491_RAW.tar", exdir="data")

#Read the CEL files into the oligo format
celpath="/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/GSE8491/data/"
list=list.files(celpath,full.names = TRUE)
data = read.celfiles(list)

#Phenodata, contains labels for the samples
ph = data@phenoData
ph
ph@data

#PARAMETERS TO CHANGE
#How to add annotation to phenoData
ph@data[ ,1] = c("NIW1","NIW2","NIW3","NIW4","NIW5","PVWKY1","PVWKY2","PVWKY3","PVWKY4","PVWKY5","PVWS1","PVWS2","PVWS3","PVWS4","PVWS5","SOW1","SOW2","SOW3","SOW4","SOW5","SOC1","SOC2","SOC3","SOC4","SOC5")  
ph
ph@data[ ,2] = c("strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","strain","control","control","control","control","control")
colnames(ph@data)[2]="source"
ph@data
groups = ph@data$source
f = factor(groups,levels=c("strain","control"))
design = model.matrix(~ 0 + f)
colnames(design) = c("strain","control")
contrast.matrix = makeContrasts(strain-control,levels=design)
Title="WKY vs SHR Juvenile - Normalised reads PCA"
Title="Normalised reads PCA"
Group="Strain"
Output="WKYVSSHR-PCA.pdf"
data.rma = rma(data)
data.matrix = exprs(data.rma)
colnames(data.matrix)=c("NIW1","NIW2","NIW3","NIW4","NIW5","PVWKY1","PVWKY2","PVWKY3","PVWKY4","PVWKY5","PVWS1","PVWS2","PVWS3","PVWS4","PVWS5","SOW1","SOW2","SOW3","SOW4","SOW5","SOC1","SOC2","SOC3","SOC4","SOC5")
m=30
mymat <- replicate(m,ph@data)

#Producing DE genes, with P<0.001 and a LOGFC 
data.fit = lmFit(data.matrix,design)

data.fit.con = contrasts.fit(data.fit,contrast.matrix)

data.fit.eb = eBayes(data.fit.con)

options(digits=2)
tab = topTable(data.fit.eb,coef=1,number=1000,adjust.method="BH")
head(tab)
topgenes = tab[tab[, "adj.P.Val"] < 0.001, ]
dim(topgenes)
topups = topgenes[topgenes[, "logFC"] > 1, ]
dim(topups)
topdowns = topgenes[topgenes[, "logFC"] < -1, ]
dim(topdowns)

DEresults = decideTests(data.fit.eb,method='global',adjust.method="BH",p.value=0.05,lfc=1)
DEresults[1:10,]

data.matrix.up = data.matrix[(rownames(topups)),]
data.matrix.down = data.matrix[(rownames(topdowns)),]

#UPREGULATED
#Creating gene symbols for the heatmap
up <- select(rae230a.db, keys = rownames(data.matrix.up), columns = c("ENTREZID", "ENSEMBL", 
                                                                      "SYMBOL"), keytype = "PROBEID")
#length(is.na(down$SYMBOL))
#length(is.na(res$SYMBOL))
#length(is.na(up$SYMBOL))

#Making the rownames into the first column
d <- data.matrix.up
names <- rownames(d)
rownames(d) <- NULL
res1 <- cbind(names,d)
colnames(res1)[colnames(res1)=="names"] <- "PROBEID"
head(res1)
res1<-data.frame(res1)

#merging the dataframes:
res2 <- data.frame(up[match(res1$PROBEID, up$PROBEID),],
                   res1)

#Removing columns:
genesymbol<-res2
genesymbol$PROBEID.2=NULL
genesymbol$PROBEID=NULL
genesymbol$ENTREZID=NULL
genesymbol$ENSEMBL=NULL
genesymbol$PROBEID.1=NULL
#genesymbol$SYMBOL=NULL

#REMOVING THE NA
matrix<-na.omit(genesymbol)

#CREATING ROWNAMES
rownames(matrix) <- NULL
rownames(matrix)=matrix$SYMBOL
matrix$SYMBOL=NULL

#Wrapping the data 
matrix2 = data.matrix(matrix)

#Heatmap: 
library(heatmaply)
heatmaply(matrix2, dendrogram = "none", Rowv = FALSE, Colv = FALSE,xlab = "samples", ylab = "genes", 
          main = "Upregulated genes",margins = c(6, 10),
          k_col = 2, k_row = 2)

#downregulated
#THIS IsN'T WORKING ATM. There are duplicates of gene symbols
down<-AnnotationDbi::select(rae230a.db, keys = rownames(data.matrix.down), columns = c("ENTREZID", "ENSEMBL", 
                                                                                       "SYMBOL"), keytype = "PROBEID")
#Making the rownames into the first column
d <- data.matrix.down
names <- rownames(d)
rownames(d) <- NULL
res1 <- cbind(names,d)
colnames(res1)[colnames(res1)=="names"] <- "PROBEID"
head(res1)
res1<-data.frame(res1)

#merging the dataframes:
res2 <- data.frame(down[match(res1$PROBEID, down$PROBEID),],
                   res1)
#Removing columns:
downgenesymbol<-res2
downgenesymbol$PROBEID.2=NULL
downgenesymbol$PROBEID=NULL
downgenesymbol$ENTREZID=NULL
downgenesymbol$ENSEMBL=NULL
downgenesymbol$PROBEID.1=NULL
#genesymbol$SYMBOL=NULL

#REMOVING THE NA
matrixdown<-na.omit(downgenesymbol)

#CREATING ROWNAMES
rownames(matrixdown) <- NULL
rownames(matrixdown)=matrixdown$SYMBOL
matrixdown$SYMBOL=NULL

#Wrapping the data 
matrixdown2 = data.matrix(matrixdown)

#Heatmap: 
library(heatmaply)
heatmaply(matrixdown2, dendrogram = "none", Rowv = FALSE, Colv = FALSE,xlab = "samples", ylab = "genes", 
          main = "Downregulated genes",margins = c(6, 10),
          k_col = 2, k_row = 2)



#Trying to use Biomart in R to obtain the gene symbols, to reduce NAs found in the data. BUT: 
#The below commands produce this error, when trying to convert probe IDs into gene symbols using biomart in R
#Error: 'with_affy_rae230a' is a boolean filter and needs a corresponding logical value of 
#TRUE or FALSE to indicate if the query should retrieve all data that fulfill the boolean or 
#alternatively that all data that not fulfill the requirement should be retrieved
require("biomaRt")
library(org.Rn.eg.db)

data.matrix.up
ens1 <- data.matrix.up
names <- rownames(ens1)
rownames(ens1) <- NULL
ens2 <- cbind(names,ens1)
colnames(ens2)[colnames(ens2)=="names"] <- "gene"
head(ens2)
ens2<-data.frame(ens2)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("rnorvegicus_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="with_affy_rae230a",
  values = list(ens2$gene, TRUE),
  uniqueRows=TRUE
  
#annotLookup1 <- data.frame(full1[match(annotLookup$ensembl_gene_id, row.names(full1)),],annotLookup)
#annotLookup1$ensembl_gene_id=NULL
#annotLookup2=annotLookup1[,c("external_gene_name","logFC","logCPM","F","PValue","FDR")]
  
#Summary of session
write.table(sessionInfo(), sep = "\t", file = "Session Info")