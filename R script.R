# load packages
#install.packages("Matrix")
#install.packages("lattice")
#install.packages("fdrtool")
#install.packages("rpart")
#install.packages("ggplot2")
#install.packages("affyPLM")
#Install Bioconductor packages
#BiocManager::available("")
#BiocManager::install(c("limma"))
#BiocManager::install(c("Biostrings"))
#BiocManager::install(c("Biobase"))
#BiocManager::install(c("genefilter"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("affyExpressionFS")
#BiocManager::install(c("'pd.rta.1.0'"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#Set working directory
setwd("/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/GSE8491/")

#Load in the packages
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)
library(ggplot2)
library(affyPLM)
library(limma)
library(Biostrings)
library(Biobase)
library(genefilter)
library(oligo)
library(ggplot2)
library(grid)
library(gridExtra)
library(DESeq2)
library(simpleaffy)
library(pd.rta.1.0)
library(rae230a.db)
library(AnnotationDbi)

# load series and platform data from GEO
library(GEOquery)
getGEOSuppFiles("GSE8491")

untar("GSE8491/GSE8491_RAW.tar", exdir="data")
#cels <- list.files("data/", pattern = "[gz]")
#length(cels)
#sapply(paste("data", cels, sep="/"), gunzip)
#cels

#Read the CEL files into the oligo format
celpath="/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/GSE8491/data/"
list=list.files(celpath,full.names = TRUE)
data = read.celfiles(list)

#Open the CEL-files in R, using affy
#Generates an AffyBatch 
#celpath ="/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/microarray"
#data = ReadAffy(celfile.path=celpath)

#Open the CEL files in R, using oligo
#Generates a HTAFeatureSet
#celpath ="/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/microarray"
#list = list.files(celpath,full.names=TRUE)
#data = read.celfiles(list)

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
res <- select(rae230a.db, keys = rownames(data.matrix), columns = c("ENTREZID", "ENSEMBL", 
                                                                    "SYMBOL"), keytype = "PROBEID")
colnames(data.matrix)=c("NIW1","NIW2","NIW3","NIW4","NIW5","PVWKY1","PVWKY2","PVWKY3","PVWKY4","PVWKY5","PVWS1","PVWS2","PVWS3","PVWS4","PVWS5","SOW1","SOW2","SOW3","SOW4","SOW5","SOC1","SOC2","SOC3","SOC4","SOC5")
m=30
mymat <- replicate(m,ph@data)

#How to list the names of the rows or columns of a matrix
rownames(data)
columnnames(data)

#Retrieving intensities using oligo 
expr = exprs(data)
int = oligo::intensity(data)

#They contain intensities of all probes so we will limit ourselves to looking at the first five rows of the resulting matrices as specified by 1:5 (the colon denotes a range in R)
expr[1:5,]
int[1:5,]
#retrieve intensities of PM probes of specific rows in the CEL files
oligo::pm(data)[1:5,]

#You can use this probe number to check if exprs() and pm() return the same values.
expr[10,]
oligo::pm(data)[2,]

#Retrieving probe annotation using affy/Oligo
feat = data@featureData
feat
feat@data

#Retrieving experiment annotation using affy/oligo
exp = data@experimentData
exp

#retrieve the name of the CDF file associated with the arrays 
#cdfName(data)

#How to retrieve the IDs of the probe sets that are represented on the arrays
#featureNames(data)

#How to retrieve the number of probe sets represented on the arrays ?
#length(featureNames(data))

#How to retrieve the number of probes represented on the arrays
#length(probeNames(data))

#On Affymetrix arrays you should see an average of 10-20 probes per gene. 

#Quality control of microarray data
#Microarray pictures
#Microarray pictures can show large inconsistencies on individual arrays.
for (i in 1:m)
{
  wJ1 = paste("image",i,".jpg",sep="")
  jpeg(wJ1)
  image(data[,i],main=mymat@data$sample[i])
  dev.off()
}

#However, if you have a set of small arrays (like the ATH arrays) and you want to plot them on a single plot, you can use the following code for the plotting:
#op = par(mfrow = c(6,6))
#for (i in 1:12){image(data[,i],main=ph@data$sample[i])}

#Creating a pseudo-chip image
Pset = fitProbeLevelModel(data)

#based on weights
#Weights represent how much the original data contribute to the model: 
#outliers are strongly downweighted because they are so different from the ideal data. 
#Weights have values between 0 and 1. 
for (i in 1:m)
{
  WEIGHTS = paste("pseudoimage",i,".jpg",sep="")
  jpeg(WEIGHTS)
  image(Pset,which=i,main=mymat$sample[i])
  dev.off()
}

#Based on residuals
#Residuals are the second quantity used for chip pseudo-images. 
# They represent the difference between the original data and the ideal data according to the model.
# So the more a residual deviates from 0, the larger the difference between the data and the model.
#Residuals can be
#positive: the intensity of the probe is larger than the ideal value according to the model
#negative: the intensity of the probe is smaller than the ideal value according to the model
#Creating pseudoimages based on residuals is slightly different from affy:
for (i in 1:m)
{
  RESIDUALS = paste("pseudoimage",i,".jpg",sep="")
  jpeg(RESIDUALS)
  image(Pset,which=i,type="pos.residuals",main=mymat$sample[i])
  dev.off()
}

#Producing Histograms
for (i in 1:m)
{
  his = paste("histogram",i,".jpg",sep="")
  jpeg(his)
  hist(target='core',data[,i],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=mymat$sample[i])
  dev.off()
}

#Creating histogram on the same page
#op = par(mfrow = c(2,6))
#for(i in 1:12){hist(target='core',data[,i],lwd=2,which='pm',type='probeset', ylab='Density',xlab='Log2 intensities',main=ph@data$sample [i])}

#creating one plot containing the histograms of all samples, wt and mutant replicates in a different color 

#color=c('green','green','green','green','green','green','red','red','red','red','red','red')
#hist(target='core',data[,1:12],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data')

#create a plot containing the histograms of all samples (each sample in a different color) using ggplot
pmexp = oligo::pm(data)
#We create two empty vectors that will serve as the two columns of the data frame:
#one to store the sample names in, called sampleNames
#one to store the log intensities in , called logs
#If you run into problems when creating these vectors and you have to repeat the plotting, always reinitialize these vectors (recreate them so that they are empty)

sampleNames = vector()
logs = vector()
for (i in 1:m)
{
  sampleNames = c(sampleNames,rep(mymat[i],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
}
logData = data.frame(logInt=logs,sampleName=sampleNames)
dataHist2 = ggplot(logData, aes(logInt, colour = sampleName)) 
dataHist2 + geom_density()


#BOXPLOTS 
#Boxplots and histograms show the same differences in probe intensity behavior between arrays
#Box plots show:
#plotting the box plots of the raw data of each array
name = "boxplot.jpg"
jpeg(name)
boxplot(target='core',data,which='pm',col='red',names=ph@data$sample) 
dev.off()

#creating the box plot of the raw data using ggplot 

sampleNames = vector()
logs = vector()
for (i in 1:m)
{
  sampleNames = c(sampleNames,rep(mymat[i],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
}
logData = data.frame(logInt=logs,sampleName=sampleNames)
dataBox = ggplot(logData,aes(sampleName,logInt))
dataBox + geom_boxplot()

#creating a boxplot of normalised intensities
#data.rma = rma(data2)
#data.matrix = exprs(data.rma)

normal= "boxplotnorm.jpg"
jpeg(normal)
boxplot(data.matrix,col='red',names=ph@data$sample)
dev.off()

#creating a box plot of normalized intensities using ggplot
#create two vectors:
#one to store the sample names in, called sampleNames
#one to store the normalized log intensities in , called normlogs
sampleNames = vector()
normlogs = vector()
for (i in 1:m)
{
  sampleNames = c(sampleNames,rep(mymat[i],dim(data.matrix)[1]))
  normlogs = c(normlogs,data.matrix[,i])
}
normData = data.frame(norm_logInt=normlogs,sampleName=sampleNames)
dataBox = ggplot(normData, aes(sampleName,norm_logInt))
dataBox + geom_boxplot() + ylim(2,16) + ggtitle("after normalization")

#Additionally, we have to make the box plot of the raw data again this time setting the Y-scale from 2 to 16.
dataBox = ggplot(logData,aes(sampleNames,logInt))
dataBox + geom_boxplot() + ylim(2,16) + ggtitle("before normalization")

#MA plots of each array

for (i in 1:m)
{
  name = paste("MAplot",i,".jpg",sep="")
  jpeg(name)
  MAplot(data,which=i)
  dev.off()
}

#MA plot of the normalised intensities
#When you compare this plot to the one created for the raw intensities, you see a much more symmetric and even spread of the data indicating that the dependence of the variability on the average expression level is not as strong as it was before normalization.
for (i in 1:m)
{
  name = paste("MAplotnorm",i,".jpg",sep="")
  jpeg(name)
  MAplot(data.rma,which=i)
  dev.off()
}

#The standard method for normalization is RMA. 
#RMA is one of the few normalization methods that only uses the PM probes
#Background correction to correct for spatial variation within individual arrays
#Log transformation to improve the distribution of the data:
#Quantile normalization to correct for variation between the arrays
#Probe normalization to correct for variation within probe sets

# RMA normalization
#eset <- oligo::rma(data2)
#data.matrix = exprs(eset)

# ... also through an alternative route
#norm.data2 <- oligo::fitProbeLevelModel(data2)

## convert oligoPLM object to ExpressionSet!
#norm.data2 <- opset2eset(norm.data2)

#PCA PLOTS
#Parameters to change

#Title="WKY vs SHR Juvenile - Normalised reads PCA"
#Title="Normalised reads PCA"
#Group="Strain"
#Output="WKYVSSHR-PCA.pdf"

#Input data

#raw_counts = read.table("shrvwkyRN6.counts", stringsAsFactors=F)
#colnames(raw_counts)=c("Gene","WKY","WKY","WKY","SHR","SHR","SHR")
#colnames(data.matrix)=c("S1","S2","S3","S4","S5","S6","W1","W2","W3","W4","W5","W6")
#head(data.matrix)
#df <- raw_counts
df1 <- data.matrix
#df <- as.data.frame(raw_counts)
df1 <-as.data.frame(data.matrix)
#row.names(df) <- paste(df$Gene, row.names(df), sep="_") 
row.names(df1) <- paste(df1$Gene, row.names(df1), sep="_") 
#df$Gene <- NULL
df1$Gene <- NULL
head(df)
#df2=t(df)
df3=t(df1)
#head(df3)
#df_pca <- prcomp(df)
df1_pca <- prcomp(df1)
#plot(df_pca$rotation[,1], df_pca$rotation[,2])
plot(df1_pca$rotation[,1], df1_pca$rotation[,2])
#df_out <- as.data.frame(df_pca$rotation)
df1_out <- as.data.frame(df1_pca$rotation)
#df_out$Group <- sapply( strsplit(as.character(row.names(df1)), "_"), "[[", 1 )
#head(df_out)
df1_out$Group <- sapply(strsplit(as.character(row.names(df3)), "_"), "[[", 1 )
#head(df1_out)

#Add the percentage of PCA
#percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- round(df1_pca$sdev / sum(df1_pca$sdev) * 100, 2)
#percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
percentage <- paste( colnames(df1_out), "(", paste( as.character(percentage), "%", ")", sep="") )

#Generate the plot
P<-ggplot(df1_out,aes(x=PC1,y=PC2, color=Group))+ xlab(percentage[1]) + ylab(percentage[2])+geom_point(size=6)+ggtitle(Title)+theme(axis.text=element_text(size=12))
axis.title.x=element_text(size=14,face="bold", margin = margin(t = 20, r = 0, b = 0, l = 0))
P     

#Save the plot as PDF
pdf(Output,width = 10,height = 10)
library(gridExtra)
yy <- grid.arrange(P,nrow=1)
op <- par(no.readonly=TRUE)
par(op)
dev.off()

#Identification of DE genes is carried out via the Limma package
#Limma uses the output of the rma() method (data.rma) as input. 
#This means that all following code is valid for all normalized Affymetrix data regardless of the package that was used for normalization
#Using Limma for comparing to groups of samples
#ph@data[ ,2] = c("strain","strain","strain","strain","strain","strain","control","control","control","control","control","control")
#colnames(ph@data)[2]="source"
#ph@data

#Now you need to tell limma which sample belongs to which group. This info is contained in the second column (the column that we called source) of the PhenoData.
#groups = ph@data$source

#The names of the groups have to be transformed into factors. 
#Factors in statistics represent the variable that you control: which sample belongs to which group. 
#You can factorize the grouping variable by using the factor() method
#f = factor(groups,levels=c("strain","control"))

#Then you need to create a design matrix, a matrix of values of the grouping variable. 
#ANOVA needs such a matrix to know which samples belong to which group. 
#Since limma performs an ANOVA or a t-test (which is just a specific case of ANOVA), it needs such a design matrix. 
#You can create it using the model.matrix() method.

#The argument of the model.matrix method is a model formula. 
#The argument of the model.matrix method is a model formula. 
#The tilde (~) in the argument specifies the right hand side of a model equation. 
#If you define the design matrix with ~0, limma will simply calculate the mean expression level in each group.

#design = model.matrix(~ 0 + f)
#colnames(design) = c("strain","control")

#Essentially, limma will compare the mean expression level in the control samples to the mean expression level in mutant samples and this for each gene. 
#So first of all, limma needs to calculate the mean expression levels using the lmFit() method. 
#This method will fit a linear model (defined in design) to the data to calculate the mean expression level in the control and in the mutant 
#The first column contains the mean log expression in control samples.
#The second column contains the mean log expression in mutant samples.

data.fit = lmFit(data.matrix,design)
data.fit$coefficients[1:10,]

#Now you have to tell limma which groups you want to compare.
#For this you define a contrast matrix defining the contrasts (comparisons) of interest by using the makeContrasts() method. 
#Using the column names of the design matrix you specify the baseline sample (control) you want to compare to and the sample of interest (mutant).

#contrast.matrix = makeContrasts(strain-control,levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)

#Now limma is ready to perform the statistical test to compare the groups. Since we are comparing two groups, it will perform a t-test.
#Limma does not perform a standard t-test but it will calculate a moderated t-statistic with shrunken standard deviation for each gene (see slides).
#An empirical Bayes method is used to shrink the variance of each gene towards a common value for all the genes.
#This is to lower the influence of very low or very high standard deviations on the t-test. 
#Since the number of replicates is very low, the standard deviations will not be very reliable, ordinary t-statistics are not recommended. 
#However, although there are only a few replications for each gene, the total number of measurements is very large. 
#So information is combined across the genes (i.e., genome-wide shrinkage) to improve performance.
#The moderated t-test is performed by using the eBayes() method.

data.fit.eb = eBayes(data.fit.con)

#The eBayes() method returns a data frame data.fit.eb containing multiple slots. To get an overview of all the slots, use the names() method:
names(data.fit.eb)
#This slot contains the coefficient of the contrast: it's simply the difference between the mean log expression in mutant samples and the mean log expression in control samples. 
#This difference is usually called a log fold change.
data.fit.eb$coefficients[1:10,]

#The eBayes() method has performed a moderated t-test on each gene. 
#That means it has calculated a t-statistic and a corresponding p-value for each gene. 
#These values can be found in the t and the p.valuje slot. 
#To view the t-statistics and p-values of the moderated t-test for the first 10 probe sets:

data.fit.eb$coefficients[1:10,]
data.fit.eb$p.value[1:10,]
data.fit.eb$t[1:10,]
#Generating a volcano plot
#A Volcano plot is generated by using the volcanoplot() method on the output of the moderated t-test.
name = "Volcano.jpg"
jpeg(name)
volcanoplot(data.fit.eb,coef=1,highlight=100)
dev.off()

#adjust for multiple testing for a single comparison
#The number of genes that is to be returned is specified by the number argument. In most cases, a ranking of genes according to evidence for DE is sufficient because only a limited number of DE genes can be used for further study.
#Additionally, the topTable() method will adjust the p-value obtained from the moderated t-test for multiple testing. The adjustment method is defined by the adjust.method argument.
#In this case, the adjustment is done using BH which is Benjamini and Hochberg's method to control the FDR.
#The meaning of the adjusted p-value is as follows: if you select all genes with adjusted p-value below 0.05 as DE, then the expected proportion of false positives in the selected group should be less than 5%. So if you select 100 DE genes at a false discovery rate of 0.05, only 5 of them will be false positives.
options(digits=2)
tab = topTable(data.fit.eb,coef=1,number=1000,adjust.method="BH")
head(tab)

#selecting a set of genes with adjusted p-values below a threshold
#You have created a subset of the table generated by topTable() with adjusted p-values (last column of tab called adj.P.Val) below a threshold (in this example the threshold is set at 0.001):
topgenes = tab[tab[, "adj.P.Val"] < 0.001, ]
dim(topgenes)
#If you want to distinguish between up- and downregulated genes and you want to include a log fold change threshold, you need to create another subset of the remaining table
topups = topgenes[topgenes[, "logFC"] > 1, ]
dim(topups)
topdowns = topgenes[topgenes[, "logFC"] < -1, ]
dim(topdowns)

#How to adjust for multiple testing for multiple comparisons 
#For each gene and each comparison it will generate the following output:
#-1: significantly downregulated
#0: no significant evidence of differential expression
#1: significantly upregulated
DEresults = decideTests(data.fit.eb,method='global',adjust.method="BH",p.value=0.05,lfc=1)
DEresults[1:10,]

#There's an alternative method for finding DE genes that allows you to specify a false discovery rate and a threshold for the log fold change in a single command.
#It is especially interesting when you do multiple comparisons, as in the example of the 3 groups of mice samples. 
#If you use topTable() for this, you have to perform the adjustment thrice, once for each comparison, each time changing the value of the coef parameter.

#How to create a list of IDs of genes that are DE according to topTable() 
IDs.up = rownames(topups)
IDs.down = rownames(topdowns)
head(IDs.down)
head(IDs.up)

#How to create a list of IDs of genes that are DE according to decideTests() ?
#ups = subset(DEresults, DEresults[,1]==1)
#downs = subset(DEresults, DEresults[,1]==-1)
#IDs.up = rownames(ups)
#IDs.down = rownames(downs)

#How to write the probe set IDs of the DE genes to your computer ?
#write.csv(IDs.up,row.names=FALSE,col.names=FALSE,quote=FALSE,file="/IDs")
#write.csv(IDs.down,row.names=FALSE,col.names=FALSE,quote=FALSE,file="")
write.csv(data.fit.eb,file = "data.fit.csv")
write.csv(IDs.down,file = "down.csv")
write.csv(res,file = "symbol.csv")

#HEATMAP:
data.matrix.up = data.matrix[(rownames(topups)),]
data.matrix.down = data.matrix[(rownames(topdowns)),]
library(heatmaply)
heatmaply(data.matrix.up, dendrogram = "none", Rowv = FALSE, Colv = FALSE,xlab = "samples", ylab = "genes", 
          main = "Upregulated genes",margins = c(6, 10),
          k_col = 2, k_row = 2)

#ALTERNATIVE HEATMAP
#Creating a heat Map: 
#data.matrix.up = data.matrix[(rownames(topups)),]
#data.matrix.down = data.matrix[(rownames(topdowns)),]
#sampleNames = vector()
#featureNames = vector()
#heatlogs = vector()
#for (i in m)
#{
#sampleNames = c(sampleNames,rep(mymat[i],dim(topdowns)[1]))
#featureNames = c(featureNames,rownames(data.matrix.down[1:dim(topdowns)[1],]))
#heatlogs = c(heatlogs,data.matrix.down[1:dim(topdowns)[1],i])
#}
#heatData = data.frame(norm_logInt=heatlogs,sampleName=sampleNames,featureName=featureNames)

#dataHeat = ggplot(heatData, aes(sampleName,featureName))
#dataHeat + geom_tile(aes(fill=norm_logInt)) + scale_fill_gradient(low="blue", high="yellow")

#Creating a Venn diagram: 
vennDiagram(DEresults)

