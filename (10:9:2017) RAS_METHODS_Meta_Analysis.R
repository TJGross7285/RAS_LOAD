####################### (10/9/2017) Generate Spearman Rho Correlation Coefficients and Additional Statistics Between All Metabolites Indexed in Amrita, Stewart, and Brian's Versions of RAS Data and Produce Related PCA Plots of Data Sets   

setwd("~/Desktop")
library(dplyr)

####Read in RAS Validation data
read_validation<-read.csv("RAS_VALIDATION_NM_1_23_transposed.csv", check.names=FALSE)

#### Generate label indicating data set of origin (RAS_Validation)
ras_label<-c("R/OCAS Validation")
DataSet<-rep(ras_label, dim(read_validation)[1])

#### Combine data set of origin, collapsed group label, and imported data frame for RAS_Validation
RAS_Validation_Final<-cbind(DataSet, read_validation)

####Read in RAS Discovery data
read_discovery<-read.csv("RAS_DISCOVERY_NM_1_23_transposed.csv", check.names=FALSE)

#### Generate label indicating data set of origin (RAS_Validation)
ras_label<-c("R/OCAS Discovery")
DataSet<-rep(ras_label, dim(read_discovery)[1])

#### Combine data set of origin, collapsed group label, and imported data frame for RAS_Discovery
RAS_Discovery_Final<-cbind(DataSet, read_discovery)

#########Merge RAS Discovery and Validation
common_cols_fpass<- intersect(colnames(RAS_Discovery_Final), colnames(RAS_Validation_Final))
Combined_Data<- rbind(subset(RAS_Validation_Final, select = common_cols_fpass), subset(RAS_Discovery_Final, select = common_cols_fpass), make.row.names=FALSE)
origin_label<-c("Amrita_Samples")
Origin<-rep(origin_label, dim(Combined_Data)[1])
Origin_Combined_Data<-cbind(Origin, Combined_Data)

##########Read in Stewart's Samples
read_stewart<-read.csv("Stewart_Data_10_5.csv", check.names=FALSE)
colnames(read_stewart)[c(1,2,3)]<-c("SampleID","DataSet","DiseaseState")
levels(read_stewart$DataSet) <- list("R/OCAS Discovery "=c("Discovery"), "R/OCAS Validation"=c("Validation"))
read_stewart$SampleID<-gsub("-","_",read_stewart$SampleID)
read_stewart$SampleID<-gsub("UCI","",read_stewart$SampleID)
read_stewart$SampleID<-gsub("RAS","",read_stewart$SampleID)
origin_label<-c("Stewart_Samples")
Origin<-rep(origin_label, dim(read_stewart)[1])
total_stewart<-cbind(Origin,read_stewart)

##########Read in Brian's Samples
brian_raw<-read.csv("Brian_Data_10_16.csv",check.names=FALSE)
Samples<-brian_raw%>%filter(`Sample Type`=="Sample")
reduced<-Samples[,c(4,11,12,19:dim(Samples)[2])]
colnames(reduced)[c(1,2,3)]<-c("SampleID","DataSet","DiseaseState")
levels(reduced$DataSet) <- list("R/OCAS Discovery "=c("Discovery"), "R/OCAS Validation"=c("Validation"))
reduced$SampleID<-gsub("-","_",reduced$SampleID)
reduced$SampleID<-gsub("UCI","",reduced$SampleID)
reduced$SampleID<-gsub("RAS","",reduced$SampleID)
origin_label<-c("Brian_Samples")
Origin<-rep(origin_label, dim(reduced)[1])
total_brian<-cbind(Origin,reduced)

###########Merge Amrita, Brian, and Stewart Samples
common_cols_origin<-intersect(colnames(Origin_Combined_Data), colnames(total_stewart))
merge<-intersect(common_cols_origin, colnames(total_brian))
All_Origin_Combined<-rbind(subset(Origin_Combined_Data, select = merge), subset(total_stewart, select = merge), subset(total_brian, select=merge), make.row.names=FALSE)
filtered_total<-dplyr::filter(All_Origin_Combined, DiseaseState!="Super")
levels(filtered_total$DiseaseState) <- list("MCIAD"=c("A"), "Converter"=c("B"), "Control"=c("normal","C"))
filtered_total$SampleID<-as.factor(filtered_total$SampleID)

####Write filtered_total to .csv file
write.csv(filtered_total,file="Stewart_Amrita_Brian_Merged_10_16.csv")

##########################
##########################

####Read-in data as generated in "check.csv" as corrected by MM in "Copy of Check 10.6.17.xlsx"
####"Stewart_Amrita_Merged_10_11.csv"%>%"Copy of Check 10.6.17.xlsx"%>%"Stewart_Amrita_Brian_Merged_10_16_Corrected.csv"
read_corrected<-read.csv("Stewart_Amrita_Brian_Merged_10_16_Corrected.csv",check.names=FALSE)
read_corrected_dropped<-read_corrected[,2:dim(read_corrected)[2]]

####Subset data into Amrita/Stewart samples and order ascending by SampleID  
Stewart<-read_corrected_dropped %>% filter(Origin=="Stewart_Samples") %>% arrange(SampleID)
Amrita<-read_corrected_dropped %>% filter(Origin=="Amrita_Samples") %>% arrange(SampleID)
Brian<-read_corrected_dropped%>%filter(Origin=="Brian_Samples")%>%arrange(SampleID)

####Check to see that list of sample IDs in two subsets of data are identical ****SHOULD BE TRUE****
identical(Stewart$SampleID,Amrita$SampleID,Brian$SampleID)

####Subset Amrita/Stewart data sets to include only metabolite abundances
amrita_metabolites<-Amrita[,5:dim(Amrita)[2]]
stewart_metabolites<-Stewart[,5:dim(Stewart)[2]]
brian_metabolites<-Brian[,5:dim(Brian)[2]]

####Compute Spearman correlation coefficients between Amrita/Stewart samples for all metabolites 
Correlation_Amrita_Stewart<- numeric(ncol(amrita_metabolites)) 
for(i in 1:ncol(amrita_metabolites)){
               Correlation_Amrita_Stewart[i] <- cor.test(amrita_metabolites[,i], stewart_metabolites[,i],method="spearman",exact=FALSE)$estimate
}

####Compute Spearman correlation coefficients between Brian/Amrita samples for all metabolites 
Correlation_Amrita_Brian<- numeric(ncol(amrita_metabolites)) 
for(i in 1:ncol(amrita_metabolites)){
               Correlation_Amrita_Brian[i] <- cor.test(amrita_metabolites[,i], brian_metabolites[,i],method="spearman",exact=FALSE)$estimate
}

####Compute Spearman correlation coefficients between Brian/Stewart samples for all metabolites 
Correlation_Brian_Stewart<- numeric(ncol(amrita_metabolites)) 
for(i in 1:ncol(amrita_metabolites)){
               Correlation_Brian_Stewart[i] <- cor.test(brian_metabolites[,i], stewart_metabolites[,i],method="spearman",exact=FALSE)$estimate
}

####Compute associated P-values for correlation coefficients in Brian/Amrita
P.Value_Amrita_Stewart<- numeric(ncol(amrita_metabolites))
for(i in 1:ncol(amrita_metabolites)){
               P.Value_Amrita_Stewart[i] <- cor.test(amrita_metabolites[,i], stewart_metabolites[,i],method="spearman",exact=FALSE)$p.value
}

####Compute associated P-values for correlation coefficients in Brian/Amrita
P.Value_Amrita_Brian<- numeric(ncol(amrita_metabolites))
for(i in 1:ncol(amrita_metabolites)){
               P.Value_Amrita_Brian[i] <- cor.test(amrita_metabolites[,i], brian_metabolites[,i],method="spearman",exact=FALSE)$p.value
}

####Compute associated P-values for correlation coefficients in Brian/Stewart
P.Value_Brian_Stewart<- numeric(ncol(amrita_metabolites))
for(i in 1:ncol(amrita_metabolites)){
               P.Value_Brian_Stewart[i] <- cor.test(brian_metabolites[,i], stewart_metabolites[,i],method="spearman",exact=FALSE)$p.value
}

####Correct P-values for multiple testing with FDR-correction
Q.Value.Brian.Stewart<-p.adjust( P.Value_Brian_Stewart,method="fdr")
Q.Value.Amrita.Brian<-p.adjust(P.Value_Amrita_Brian,method="fdr")
Q.Value.Amrita.Stewart<-p.adjust(P.Value_Amrita_Stewart,method="fdr")

####Calculate median values for each metabolite in Amrita/Stewart data
median_amrita<-apply(amrita_metabolites,2,median,na.rm=TRUE)
median_stewart<-apply(stewart_metabolites,2,median,na.rm=TRUE)
median_brian<-apply(brian_metabolites,2,median,na.rm=TRUE)

####Calculate interquartile range for each metabolite in Amrita/Stewart data
iqr_amrita<-apply(amrita_metabolites,2,IQR,na.rm=TRUE)
iqr_stewart<-apply(stewart_metabolites,2,IQR,na.rm=TRUE)
iqr_brian<-apply(brian_metabolites,2,IQR,na.rm=TRUE)


####Generate final table with all calculated statistics 
Final_Table<-cbind(median_amrita, iqr_amrita, median_stewart, iqr_stewart, median_brian, iqr_brian, Correlation_Amrita_Stewart, Q.Value.Amrita.Stewart, Correlation_Amrita_Brian,Q.Value.Amrita.Brian,Correlation_Brian_Stewart, Q.Value.Brian.Stewart)
colnames(Final_Table)<-c("Median: Amrita Samples [uM]", "IQR: Amrita Samples [uM]","Median: Stewart Samples [uM]", "IQR: Stewart Samples [uM]", "Median: Brian Samples [uM]", "IQR: Brian Samples [uM]", "Spearman Correlation: Stewart-Amrita","Q-Value: Stewart-Amrita","Spearman Correlation: Amrita-Brian","Q-Value: Amrita-Brian", "Spearman Correlation: Stewart-Brian","Q-Value: Stewart-Brian")

write.csv(Final_Table,file="Amrita_Stewart_Correlations_10_17_Final.csv")


###########################################
###########################################

############################## (10/15/2017) PCA Plot of All Versions of RAS Data with ADRC and Supernormal Controls 

setwd("~/Desktop")
###################Read in All Versions of RAS Data
read_corrected<-read.csv("Stewart_Amrita_Brian_Merged_10_16_Corrected.csv", check.names=FALSE)
read_corrected_dropped<-read_corrected[,2:dim(read_corrected)[2]]
filter_correct<-dplyr::filter(read_corrected_dropped,DiseaseState=="Control"|DiseaseState=="Converter")
SampleID<-filter_correct$SampleID
CollapsedGroups<-filter_correct$DiseaseState
levels(CollapsedGroups) <- list("Control"=c("Control"), "Case"=c("Converter"))
abundances_RAS<-filter_correct[,5:dim(filter_correct)[2]]
DataSet<-filter_correct$Origin
levels(DataSet)<-list("RAS: Amrita"=c("Amrita_Samples"), "RAS: Brian"=c("Brian_Samples"), "RAS: Stewart"=c("Stewart_Samples"))
RAS_Final<-cbind(SampleID, DataSet, CollapsedGroups, abundances_RAS)

setwd("~/Desktop/na_replaced")
###################Read in Supernormal Controls metabolomics data and prepare for merging
read_Super<-read.csv("SUPER_NM_NA_transposed.csv",check.names=FALSE)
filter_Super<-dplyr::filter(read_Super, DiseaseState=="normal")
SampleID<-filter_Super$SampleID
CollapsedGroups<-filter_Super$DiseaseState
levels(CollapsedGroups) <- list("Control"=c("normal"))
abundances_Super<-filter_Super[,4:dim(read_Super)[2]]
Super_label<-c("Supernormal")
DataSet<-rep(Super_label, dim(filter_Super)[1])
Super_Final<-cbind(SampleID, DataSet, CollapsedGroups, abundances_Super)

###################Read in UCI ADRC metabolomics data and prepare for merging
read_ADR<-read.csv("ADR_NM_NA_transposed.csv", check.names=FALSE)
SampleID<-read_ADR$SampleID
CollapsedGroups<-read_ADR$DiseaseState
levels(CollapsedGroups) <- list("Control"=c("E"), "Case"=c("A", "B","C","D"))
abundances_ADR<-read_ADR[,4:dim(read_ADR)[2]]
ras_label<-c("UCI ADRC")
DataSet<-rep(ras_label, dim(read_ADR)[1])
ADRC_Final<-cbind(SampleID, DataSet, CollapsedGroups, abundances_ADR)


#########Merge UCI ADRC, RAS Validation, Supernormal Controls, and RAS Discovery data sets
common_cols_fpass<- intersect(colnames(ADRC_Final), colnames(RAS_Final))
common_cols_final<-intersect(common_cols_fpass,colnames(Super_Final))
Combined_Data<- rbind(subset(RAS_Final, select = common_cols_final), subset(ADRC_Final, select = common_cols_final), subset(Super_Final, select=common_cols_final), make.row.names=FALSE)

#####Extract metabolite abundances and log transform (replace INF values with NA) 
pca.abunds<-Combined_Data[,4:dim(Combined_Data)[2]]
logged<-log2(pca.abunds)
inf.index<-apply(logged,2,is.infinite)
logged[inf.index]<-NA

#####Drop metabolites with missingness in excess of 33.3%
nas<-apply(logged,2,is.na)
sum.nas<-apply(nas,2,sum)
index<-sum.nas/dim(logged)[1]*100>33.3
Filtered_Combined_Data<-logged[,index==FALSE]

#####Remove Outliers and print details regarding outlier samples
outliers<-c(76,28,259,286,13,49,67,178,362)
write.csv(Filtered_Combined_Data[outliers,],file="Outlier_RAS_All_Versions_Super_ADRC.csv")
outliers_removed<-Filtered_Combined_Data[-outliers,]

#####Perform PCA on metabolite abundances 
library(mixOmics)
pc<-pca(outliers_removed,ncomp=2)
plotIndiv(pc,comp=c(1,2),group=Combined_Data$DataSet[-outliers],legend=TRUE,ind.names=Combined_Data$SampleID[-outliers])

#####Print PDF of PCA plot
pdf(file="RAS_All_Versions_Super_ADRC_PCA.pdf")
plotIndiv(pc,comp=c(1,2),group=Combined_Data$DataSet[-outliers],legend=TRUE,ind.names=Combined_Data$SampleID[-outliers])
dev.off()

#####Create data table of metabolite loadings on principle components
loadings<-data.frame(matrix(unlist(pc$loadings), nrow=length(colnames(outliers_removed)), byrow=F))
final_loadings<-cbind(colnames(outliers_removed),as.data.frame(loadings))
colnames(final_loadings)<-c("Metabolite","PC1","PC2")

#####Write data table to csv file
write.csv(final_loadings, file="PCA_Loadings_RAS_Super_ADRC.csv")

############################## (10/15/2017) PCA Plot of All Versions of RAS Data 
setwd("~/Desktop")
#####Read in All Versions of RAS Data
read_corrected<-read.csv("Stewart_Amrita_Brian_Merged_10_16_Corrected.csv", check.names=FALSE)
read_corrected_dropped<-read_corrected[,2:dim(read_corrected)[2]]
filter_correct<-dplyr::filter(read_corrected_dropped,DiseaseState=="Control"|DiseaseState=="Converter")
SampleID<-filter_correct$SampleID
CollapsedGroups<-filter_correct$DiseaseState
levels(CollapsedGroups) <- list("Control"=c("Control"), "Case"=c("Converter"))
abundances_RAS<-filter_correct[,5:dim(filter_correct)[2]]
DataSet<-filter_correct$Origin
levels(DataSet)<-list("RAS: Amrita"=c("Amrita_Samples"), "RAS: Brian"=c("Brian_Samples"), "RAS: Stewart"=c("Stewart_Samples"))
RAS_Final<-cbind(SampleID, DataSet, CollapsedGroups, abundances_RAS)

#####Extract metabolite abundances and log transform (replace INF values with NA) 
pca.abunds<-RAS_Final[,4:dim(RAS_Final)[2]]
logged<-log2(pca.abunds)
inf.index<-apply(logged,2,is.infinite)
logged[inf.index]<-NA

#####Drop metabolites with missingness in excess of 33.3%
nas<-apply(logged,2,is.na)
sum.nas<-apply(nas,2,sum)
index<-sum.nas/dim(logged)[1]*100>33.3
Filtered_Combined_Data<-logged[,index==FALSE]

#####Remove Outliers and print details regarding outlier samples
outliers<-c(76,28,13,286,259,49,67,178)
write.csv(Filtered_Combined_Data[outliers,],file="Outlier_RAS_All_Versions.csv")
outliers_removed<-Filtered_Combined_Data[-outliers,]

#####Perform PCA on metabolite abundances 
library(mixOmics)
pc<-pca(outliers_removed,ncomp=2)
plotIndiv(pc,comp=c(1,2),group=RAS_Final$DataSet[-outliers],legend=TRUE,ind.names=RAS_Final$SampleID[-outliers])

#####Print PDF of PCA plot
pdf(file="RAS_All_Versions_PCA.pdf")
plotIndiv(pc,comp=c(1,2),group=RAS_Final$DataSet[-outliers],legend=TRUE,ind.names=RAS_Final$SampleID[-outliers])
dev.off()

#####Create data table of metabolite loadings on principle components
loadings<-data.frame(matrix(unlist(pc$loadings), nrow=length(colnames(outliers_removed)), byrow=F))
final_loadings<-cbind(colnames(outliers_removed),as.data.frame(loadings))
colnames(final_loadings)<-c("Metabolite","PC1","PC2")

#####Write data table to csv file
write.csv(final_loadings, file="PCA_Loadings_RAS_All_Versions.csv")