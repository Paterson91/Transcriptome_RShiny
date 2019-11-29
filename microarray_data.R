# load packages
install.packages("Matrix")
install.packages("lattice")
install.packages("fdrtool")
install.packages("rpart")
install.packages("ggplot2")
install.packages("affyPLM")
#Install Bioconductor packages
BiocManager::available("")
BiocManager::install(c("limma"))
BiocManager::install(c("Biostrings"))
BiocManager::install(c("Biobase"))
BiocManager::install(c("genefilter"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("simpleaffy")
BiocManager::install(c("rta10cdf"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



#Using a method from a specific package
oligo::intensity(data)

#Documentation in R studio
help(package="affy”)

#How to check which data type a variable belongs to
class(objectname)

#How to list the elements of a data frame
names(dataframename)

#How to list the names of the rows or columns of a matrix
rownames(matrixname)
columnnames(matrixname)

#Open the CEL-files in R, using oligo
# specify the path on your computer where the folder that contains the CEL-files is located
celpath ="/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/microarray"

# import CEL files containing raw probe-level data into an R AffyBatch object
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)
#Came up with error that could not find function "read.celfiles", so use these commands to correct error
library(oligo)
celFiles <- list.celfiles("/Users/bethalexander/Desktop/university/year two/Summer placement/working directory/microarray")

#Retrieving intensities using oligo 
#Need to specify oligo:: in front of the intensity() and the pm() method
expr = exprs(data)
int = intensity(data)

#They contain intensities of all probes so we will limit ourselves to looking at the first five rows of the resulting matrices as specified by 1:5 (the colon denotes a range in R)
expr[1:5,]
int[1:5,]
#retrieve intensities of PM probes of specific rows in the CEL files
pm(data)[1:5,]

#You can use this probe number to check if exprs() and pm() return the same values.
expr[10,]
pm(data)[2,]
#Now the results of the pm() and the exprs() method do correspond.
#Warning: The first column of pm(data) contains probe numbers not probe set IDs !

#How to retrieve intensities of the PM probes of a specific probe set 

pm(data)[9,10,11,12]


#plot the intensities of the PM probes of a specific probe set
#Firstly you need the data in the correct format 
#the variables you want to plot on each axis have to be in different columns.
#You want to plot the intensities (Y-axis) of each probe (probe ID on the X-axis) of a probe set. 
#So you need a matrix with probe IDs in the first column and intensities in the second column.
#In this way we will create the vector (column) containing the probe IDs.

probeNrs = rep(rownames(pm(data)[11,]),12)

#Then, we create the vector (column) containing the intensities of the probes of the probe set. To create a vector you can use the c() command.

ints=c(pm(data)[11])[1],pm(data)[11])[2],pm((data)[11])
[3],pm((data)[11])[4],pm((data)[11])[5],pm((data)[11])[6]),
(pm(data)[11])[7],(pm(data)[11])[8],(pm(data)[11])[9],(pm(data)[11])[10],(pm(data)[11])[11],(pm(data)[11])[12]

ints=c(pm(data)[11,])[,1],pm(data)[11,])[,2],pm(data,"245027_at")
[,3],pm(data,"245027_at")[,4],pm(data,"245027_at")[,5],pm(data,"245027_at")[,6]),
pm(data,"245027_at")[,1],pm(data,"245027_at")[,2],pm(data,"245027_at")
[,3],pm(data,"245027_at")[,4],pm(data,"245027_at")[,5],pm(data,"245027_at")[,6])


#Now we will combine the vector with the probe IDs and the vector with the intensities into a data frame (= matrix in which columns contain data of different types: numbers, text, boolean...). The data frame contains two columns: one is called probeNrs and one is called ints.
pset = data.frame(probeNrs=probeNrs)

#Since probe ID is a categorical variable (it can take only a limited number of different values) you have to transform it into a factor so that R knows that they are categorical and will treat them likewise. Transform probe IDs into factors
pset$PNs = factor(pset$probeNrs,levels=pset$probeNrs)

#Now the data is in the correct format for ggplot. 
#The first argument is the data you want to plot, 
#the second defines the X- and the Y-axis. 
#In the second line of code you define the symbol you want to use on the plot (geom_point()) and the titles of the X- and Y-axis.
scatter = ggplot(pset,aes(PNs,ints))
scatter + geom_point() + labs(x="Probe Number",y="Intensity of PM probe")

#plot the intensities of the PM probes of a specific probe set, 
#coloring them according to the group the sample belongs to (wild type or mutant) 

arrays = c(rep("wt",33),rep("mut",33))

#We add this vector as a column called arrays to the data frame
pset$arrays = arrays

#We can now colour the symbols according to their value in the arrays column by adding the argument colour = arrays to the ggplot() command. To create a legend for the colouring you have to add the same argument to the labs() command.
scatter = ggplot(pset,aes(PNs,ints,colour=arrays))
scatter + geom_point() + labs(x="Probe Number",y="Intensity of PM probe",colour=arrays)

#How to retrieve the sample annotation of the data
#Phenodata, contains labels for the samples

ph = data@phenoData
ph

#To look at all the data in the data frame ask for the data slot.

ph@data

#Retrieving probe annotation using affy/Oligo
#Microarray data sets should include information on the probes. 
#AffyBatches have a slot called featureData, a data frame that contains labels for the probes.
#For most data sets (also public data coming from GEO or ArrayExpress) the featureData has not been defined.

feat = data@featureData
feat
feat@data

#Retrieving experiment annotation using affy/oligo
#Microarray data sets should also include information on the experiment. 
#AffyBatches have a slot for this called experimentData. For most data sets (also public data coming from GEO or ArrayExpress) the featureData has not been defined.

exp = data@experimentData
exp
#retrieve the name of the CDF file associated with the arrays 
cdfName(data)

#How to retrieve the IDs of the probe sets that are represented on the arrays
featureNames(data)

#How to retrieve the number of probe sets represented on the arrays ?
length(featureNames(data))

#How to retrieve the number of probes represented on the arrays
length(probeNames(data))

#On Affymetrix arrays you should see an average of 10-20 probes per gene. 

#Quality control of microarray data
#How to add annotation to phenoData
ph@data[ ,1] = c("wJ1","wJ2","wJ3","wJ4","wJ4","wJ6","SJ1","SJ2","SJ3","SJ4","SJ5","SJ6")
ph

#Creating plots to assess the quality of the data
#Instead of printing these plots in RStudio or the R editor, we will save the plots to our hard drive
#Microarray pictures
#Microarray pictures can show large inconsistencies on individual arrays.
for (i in 1:12)
{
  wJ1 = paste("image",i,".jpg",sep="wJ1")
  jpeg(wJ1)
  image(data[,i],main=ph@data$sample[i])
  dev.off()
}

#However, if you have a set of small arrays (like the ATH arrays) and you want to plot them on a single plot, you can use the following code for the plotting:
op = par(mfrow = c(6,6))
for (i in 1:12){image(data[,i],main=ph@data$sample[i])}
     
#Creating a pseudo-chip image
#useful for detecting spatial differences (artifacts) on the invidual arrays (so not for comparing between arrays).
#Pseudo-images are generated by fitting a probe-level model (PLM) to the data that assumes that all probes of a probe set behave the same in the different samples: 
#probes that bind well to their target should do so on all arrays, probes that bind with low affinity should do so on all arrays.    
# fitting a probe-level model to the data
Pset = fitProbeLevelModel(data)

#based on weights
#Weights represent how much the original data contribute to the model: 
#outliers are strongly downweighted because they are so different from the ideal data. 
#Weights have values between 0 and 1. 
#So the smaller the weight of a probe
#-> the more the probe is not showing the typical behavior that it shows on the other arrays
#-> the more its intensity can be considered an outlier
#Creating pseudoimages based on weights is done exactly the same as in affy.
for (i in 1:12)
{
  WEIGHTS = paste("pseudoimage",i,".jpg",sep="")
  jpeg(WEIGHTS)
  image(Pset,which=i,main=ph@data$sample[i])
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
#The type argument is used to control which of these images is drawn:
#type="resids" gives a pseudo-image of residuals
#type="pos.resids" only high positive residuals are drawn in red, while negative and near 0 residuals being drawn in white
#type="neg.resids" only extreme negative residuals are drawn in blue, while positive negative and near 0 residuals being drawn in white
#type="sign.resids" gives images where all negative residuals regardless of magnitude are indicated by blue and all positive residuals by red

for (i in 1:12)
{
  RESIDUALS = paste("pseudoimage",i,".jpg",sep="")
  jpeg(RESIDUALS)
  image(Pset,which=i,type="residuals",main=ph@data$sample[i])
  dev.off()
}

#Producing Histograms
#Plotting histograms of each array
#The which argument allows you to specify if you want to use:
#perfect match probes only: which='pm'
#mismatch probes only: which='mm'
#both: which='both'
#creating histograms on different pages

for (i in 1:12)
{
  his = paste("histogram",i,".jpg",sep="")
  jpeg(his)
  hist(target='core',data[,i],lwd=2,which='pm',ylab='Density',xlab='Log2 intensities',main=ph@data$sample[i])
  dev.off()
  }

#Creating histogram on the same page

op = par(mfrow = c(2,6))
for(i in 1:12){hist(target='core',data[,i],lwd=2,which='pm',type='probeset', ylab='Density',xlab='Log2 intensities',main=ph@data$sample [i])}

#creating one plot containing the histograms of all samples, wt and mutant replicates in a different color 

color=c('green','green','green','green','green','green','red','red','red','red','red','red')
hist(target='core',data[,1:12],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data')

#create a plot containing the histograms of all samples (each sample in a different color) using ggplot
pmexp = pm(data)
#We create two empty vectors that will serve as the two columns of the data frame:
#one to store the sample names in, called sampleNames
#one to store the log intensities in , called logs
#If you run into problems when creating these vectors and you have to repeat the plotting, always reinitialize these vectors (recreate them so that they are empty)

sampleNames = vector()
logs = vector()
  for (i in 1:12)
  {
  sampleNames = c(sampleNames,rep(ph@data[i,1],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
  }
logData = data.frame(logInt=logs,sampleName=sampleNames)
dataHist2 = ggplot(logData, aes(logInt, colour = sampleName)) 
  dataHist2 + geom_density()

#BOXPLOTS 
#Boxplots and histograms show the same differences in probe intensity behavior between arrays
#Box plots show:
#the median: center value, half of the intensities are lower than this value, half of the intensities are higher (= line in the box)
#the upper quartile: a quarter of the values are higher than this quartile (= upper border of the box)
#the lower quartile: a quarter of the values are lower than this quartile (= lower border of the box)
#the range: minimum and maximum value (= borders of the whiskers)
#individual extreme values (= points outside the whiskers)

#plotting the box plots of the raw data of each array
name = "boxplot.jpg"
jpeg(name)
boxplot(target='core',data,which='pm',col='red',names=ph@data$sample) 
dev.off()

#creating the box plot of the raw data using ggplot 
pmexp = pm(data)
sampleNames = vector()
logs = vector()
for (i in 1:12)
{
sampleNames = c(sampleNames,rep(ph@data[i,1],dim(pmexp)[1]))
logs = c(logs,log2(pmexp[,i]))
}
logData = data.frame(logInt=logs,sampleName=sampleNames)
dataBox = ggplot(logData,aes(sampleName,logInt))
dataBox + geom_boxplot()


