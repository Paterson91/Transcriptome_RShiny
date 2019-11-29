##GeneOntology and GSEA visualization##

#Setting working directory
setwd("/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/micrarray_data/")

#Installing the packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("AnnotationHub")
#BiocManager::install("org.Rn.eg.db")
#BiocManager::install("DOSE")

#OBTAINING GENE SYMBOL FROM ENSEMBL ID
require("biomaRt")
library(org.Rn.eg.db)

#Input data
fold <- read.csv("/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/micrarray_data/fold.csv")
colnames(fold) = c("ensembl_gene_id","LogFC")

#Entrezgene
mart <-useMart("ENSEMBL_MART_ENSEMBL")
mart <-useDataset("rnorvegicus_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id","entrezgene_id"),
  values=fold$ensembl_gene_id,
  uniqueRows=TRUE)
annotLookup1 <- data.frame(fold[match(annotLookup$ensembl_gene_id, fold$ensembl_gene_id),],
                           annotLookup)

annotLookup2<-na.omit(annotLookup1)
annotLookup2$ensembl_gene_id=NULL
annotLookup2$ensembl_gene_id.1=NULL

#changing the position of the columns
annotLookup3<-annotLookup2[,c(2,1)]

#Input data
#How to make a seperate subset of the data
#m<-topgenes[,1:2]
#m$AveExpr=NULL
genelist<-annotLookup3
## feature 1: numeric vector
geneList <- annotLookup3[,2]
## feature 2: named vector
names(geneList) <- as.character(annotLookup3[,1])
## feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)
#define fold change greater than 2 as DEGs
gene <- names(geneList)[abs(geneList) > 2]
head(gene)

#GO analysis - can't get my own dataset into the right format

#Load the sample data into R
#data(geneList, package="DOSE")
#head(geneList)

#WIKIPathways analysis
#library(magrittr)
#library(clusterProfiler)
#data(geneList, package="DOSE")
#gene <- names(geneList)[abs(geneList) > 2]
#wpgmtfile <- system.file("/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/micrarray_data/rat/wikipathways-20190810-gmt-Rattus_norvegicus.gmt", package="clusterProfiler")
#wp2gene <- read.gmt(wpgmtfile)
#wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
#wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
#wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

#wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Rattus_norvegicus.gmt", package="clusterProfiler")
#wp2gene <- read.gmt(wpgmtfile)
#wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
#wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
#wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

#GO analysis
#library(clusterProfiler)
#data(geneList, package="DOSE")
#gene <- names(genelist)[abs(genelist) > 2]
#gene.df <- bitr(gene, fromType = "ENTREZID",
#             toType = c("ENSEMBL", "SYMBOL"),
#            OrgDb = org.Rn.eg.db)
#Bar plot
#library(DOSE)
#data(geneList)
#de <- names(geneList)[abs(geneList) > 2]
#edo <- enrichDGN(de)
#library(enrichplot)
#barplot(edo, showCategory=20)

#Dot plot
#edo2 <- gseNCG(geneList, nPerm=10000)
#p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
#p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
#plot_grid(p1, p2, ncol=2)

#ridgeline plot for expression distribution of GSEA result
#ridgeplot(edo2)

## convert gene ID to Symbol
#edox <- setReadable(edo, 'org.Rn.eg.db', 'ENTREZID')
#cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
#cnetplot(edox, categorySize="pvalue", foldChange=geneList)
#make it colourful
#cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

#Heatmap
#heatplot(edox)
#heatplot(edox, foldChange=geneList)

#Enrichment map
#emapplot(edo)

#UpSet plot
#upsetplot(edo)

#Running score and preranked list of GSEA result
#gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
#gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1])
#gseaplot(edo2, geneSetID = 1, title = edo2$Description[1])

#Another method to plot GSEA result is the gseaplot2 function
#gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])

#multile gene sets displayed on the same figure
#gseaplot2(edo2, geneSetID = 1:3)

#displaying the pvalue table on the plot via pvalue_table parameter:
#gseaplot2(edo2, geneSetID = 1:3, pvalue_table = TRUE,
#          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

#User can specify subplots to only display a subset of plots:
#gseaplot2(edo2, geneSetID = 1:3, subplots = 1)
#gseaplot2(edo2, geneSetID = 1:3, subplots = 1:2)

#The gsearank function plot the ranked list of genes belong to the specific gene set
#gsearank(edo2, 1, title = edo2[1, "Description"])

#Multiple gene sets can be aligned using cowplot
#library(ggplot2)
#library(cowplot)

#pp <- lapply(1:3, function(i) {
#  anno <- edo2[i, c("NES", "pvalue", "p.adjust")]
#  lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

#  gsearank(edo2, i, edo2[i, 2]) + xlab(NULL) +ylab(NULL) +
#    annotate("text", 0, edo2[i, "enrichmentScore"] * .9, label = lab, hjust=0, vjust=0)
#})
#plot_grid(plotlist=pp, ncol=1)

#pubmed trend of enriched terms
#terms <- edo$Description[1:3]
#p <- pmcplot(terms, 2010:2017)
#p2 <- pmcplot(terms, 2010:2017, proportion=FALSE)
#plot_grid(p, p2, ncol=2)

#goplot
#goplot(ego)

#To view the KEGG pathway, user can use browseKEGG function, which will open web browser and highlight enriched genes.
#browseKEGG(kk, 'hsa04110')


###GO Analysis

##ALEX'S CLUSTERPROFILE COMMANDS##

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("org.Mm.eg.db", version = "3.8")
BiocManager::install("biomaRt", version = "3.8")

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My Analysis/DESeq2 - Day vs Night/")
library(clusterProfiler)
library(enrichplot)
require("biomaRt")
library(org.Mm.eg.db)

#LOAD GENES AND ANNOTATE
rnorvegicus_gene_ensembl
mmusculus_gene_ensembl
n=11

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("rnorvegicus_gene_ensembl", mart)
ens = read.csv("Full_table.csv",header = T, row.names = 1)

if (exists("annotLookup_table")) {
print("Lookup Table already created")
} else {
annotLookup_table <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name","entrezgene_id"),
  filter="ensembl_gene_id",
  values=ens$gene,
  uniqueRows=TRUE)
annotLookup_match = merge(x=ens, y=annotLookup_table, by.x="gene", by.y = "ensembl_gene_id", all.x = TRUE)
}

head(ens)
head(annotLookup_table)
head(annotLookup_match)

GeneList_reads = annotLookup_match[,c(1,8:ncol(annotLookup_match))]
GeneList_reads$gene_biotype = NULL
GeneList_reads$external_gene_name=NULL
ordering <- order(rowMeans(GeneList_reads),decreasing=T)[1:n]
GeneList_reads_ordered = GeneList_reads[ordering,]
head(GeneList_reads_ordered)


GeneList_fc = annotLookup_match[,c(1,3)]
GeneList_fc$log2FoldChange = 2^GeneList_fc$log2FoldChange
colnames(GeneList_fc)=c("Gene","FoldChange")
GeneList_fc_order = GeneList_fc[order(-GeneList_fc$FoldChange),]
head(GeneList_fc_order)


### GO Classification

#Set levels
#Set choice of; CC - Cellular Component, MF - Molecular Function, BP - Biological Process

ont_selection = "CC"
level_selection = 7

ggo1 <- groupGO(gene     = as.character(annotLookup_match$entrezgene),
               keyType = "ENTREZID",
               OrgDb    = org.Mm.eg.db,
               ont      = ont_selection,
               level    = level_selection,
               readable = TRUE)
barplot(ggo1, drop=TRUE, showCategory=12, border = c(10,10))
dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_Go.png"))
dev.off()


### GO over-representation test


ego <- enrichGO(gene          = as.character(annotLookup_match$ensembl_gene_id),
                OrgDb         = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont           = ont_selection,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
barplot(ego, showCategory=8)
dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_OverRep.png"))
dev.off()

### GSEA test

geneList = as.numeric(GeneList_fc_order[,2])
names(geneList) = as.character(GeneList_fc_order[,1])

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "CC",
              keyType ="ENSEMBL",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

ridgeplot(ego3)

gseaplot(ego3, geneSetID = 2, title = ego3$Description[2])
gseaplot2(ego3, geneSetID = 1:6)


### Pubmed Trend of enriched terms

terms <- ego3$Description[1:3]
p <- pmcplot(terms, 2010:2017)
p2 <- pmcplot(terms, 2010:2017, proportion=FALSE)
plot_grid(p, p2, ncol=2)


# Visualisation

dotplot(ggo)
emapplot(ego)
goplot(ego)
gseaplot(kk2, geneSetID = "hsa04145")
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=geneList)


write.table(sessionInfo(), sep = "\t", file = "Session Info")
