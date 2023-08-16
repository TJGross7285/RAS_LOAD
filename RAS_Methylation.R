library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)

dataDirectory<-c("/data/users/tjgross/RAS_Methyl")
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet_10_12_2019_TJG.csv")
rgSet <- read.metharray.exp(targets=targets)


sampleNames(rgSet) <- targets$Sample_Name
detP <- detectionP(rgSet)
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
targets <- targets[keep,]
detP <- detP[,keep]

mSetSq <- preprocessQuantile(rgSet) 

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                        c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
xReactiveProbes <- read.csv("48639-non-specific-probes-Illumina450k.csv",sep="/", stringsAsFactors=FALSE)

keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,] 

mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)






use<-cbind(as.data.frame(rownames(ann450k)),as.data.frame(ann450k[,24:25]))










load("Methyl_annot_2_25.RData",verbose=TRUE)
Y<-tidyr::separate_rows(use,UCSC_RefGene_Name,sep=";")

R<-read.csv("Human Metabolic Pathways Panel Gene List.csv",check.names=FALSE)
colnames(R)[1]<-"UCSC_RefGene_Name"
join<-dplyr::inner_join(Y,R,by="UCSC_RefGene_Name")
index<-duplicated(join$Probe)
join<-join[index==FALSE,]


load("mVals_2_24_2020.RData",verbose=TRUE)
index<-rownames(mVals) %in% join$Probe 
mVals_subset<-mVals[index==TRUE,]
mVals_subset<-cbind(as.data.frame(rownames(mVals_subset)),as.data.frame(mVals_subset))
colnames(mVals_subset)[1]<-"Probe"

fjoin<-dplyr::inner_join(join,mVals_subset,by="Probe")




