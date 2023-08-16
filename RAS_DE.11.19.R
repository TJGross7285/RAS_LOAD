library(sva)
library(dplyr)
load("RAS_By_Group.RData",verbose=TRUE)
ras_pos<-rbind(converter_POS,control_POS,mciad_POS,super_POS)
ras_neg<-rbind(converter_NEG,control_NEG,mciad_NEG,super_NEG)
identical(ras_pos$SampleID,ras_neg$SampleID)

####Read in abundance data and subset to final table of features 
ras_pos<-rbind(converter_POS,control_POS,mciad_POS,super_POS)[,-c(1:2)]
ras_pos_mz<-gsub("_.+$","",colnames(ras_pos))
ras_pos_rt<-gsub("^.+_","",colnames(ras_pos))
ras_neg<-rbind(converter_NEG,control_NEG,mciad_NEG,super_NEG)[,-c(1:2)]
ras_neg_mz<-gsub("_.+$","",colnames(ras_neg))
ras_neg_rt<-gsub("^.+_","",colnames(ras_neg))


meta_pheno<-rbind(converter_POS,control_POS,mciad_POS,super_POS)[,c(1:2)]
meta_abunds<-cbind(rbind(ras_pos_mz,ras_pos_rt),rbind(ras_neg_mz,ras_neg_rt))
abunds<-log2(cbind(ras_pos,ras_neg))
mode<-as.factor(c(rep("ESI+",dim(ras_pos)[2]),rep("ESI-",dim(ras_neg)[2])))
na_index<-apply(abunds,2,is.na)
infinite_index<-apply(abunds,2,is.infinite)
final_index<-na_index==TRUE|infinite_index==TRUE
abunds[final_index==TRUE]<-NA
index<-caret::nearZeroVar(abunds)
meta_abunds<-meta_abunds[,-index]
mode<-mode[-index]
abunds<-abunds[,-index]

####Impute missing data with KNN imputation
library(impute)
T_abunds<-t(abunds)
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
edata<-apply(edata,2,as.numeric)


####Conduct SVA over samples 
pheno_table<-meta_pheno[,1:2]
colnames(pheno_table)<-c("StudyID","DiseaseState")
mod<-model.matrix(~as.factor(DiseaseState),data=pheno_table)
mod0<-model.matrix(~1,data=pheno_table)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)$sv
total<-cbind(as.data.frame(pheno_table),as.data.frame(svobj))
colnames(total)<-c("SampleID",
					"Main",
					"SV1",
					"SV2",
					"SV3",
					"SV4")
####Set up accessory objects for DE incorporating pre-post timepoints
design1<-model.matrix(~0+Main+SV1+SV2+SV3+SV4,data=total)
colnames(design1)[1:4]<-c("Control","Converter","MCIAD","Super")

####Fit linear model for pairwise contrasts 
fit1<-lmFit(edata,design1)
cm <- makeContrasts(`Converter-Control`= Converter-Control,
					`MCIAD-Converter`= MCIAD-Converter,
					`Super-Control`= Super-Control,
					`Super-MCIAD`= Super-MCIAD,
					`Super-Converter`= Super-Converter,
					`MCIAD-Control`= MCIAD-Control,  
					levels=design1)

fit1_F <- contrasts.fit(fit1, cm)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
final_DE_matrix<-cbind(rownames(T),T)
colnames(final_DE_matrix)[1]<-"Feature"

write.table(FINAL%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","MCIAD.Control"),
			sep="\t",file="RAS_POS_10_30.txt",row.names=FALSE)
write.table(FINAL%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","MCIAD.Control"),
			sep="\t",file="RAS_NEG_10_30.txt",row.names=FALSE)

cd "/Users/TGross/Desktop/Mummichog for ADRC Paper" 
source activate py27
mummichog -f RAS_POS_10_30.txt -o RAS_POS_10_30 -m positive -u 7 
mummichog -f RAS_NEG_10_30.txt -o RAS_NEG_10_30 -m negative -u 7 

write.csv(T, file="PACE_DE_Collapsed_6_25.csv")






















####Carry out SVA
edata<-imputed_log$data
mod<-model.matrix(~as.factor(factor(DiseaseState)), data=meta_pheno)
mod0<-model.matrix(~1,data=meta_pheno)
n.sv<-num.sv(edata,mod,seed=122)
svobj<-sva(edata,mod,mod0,n.sv=n.sv)
modSv<-cbind(mod,svobj$sv)
mod0Sv<-cbind(mod0,svobj$sv)

median_case<-apply(abunds[meta_pheno$DiseaseState=="Converter",],2,median)
median_control<-apply(abunds[meta_pheno$DiseaseState=="Control",],2,median)
log2FC<-log2(exp(median_case)/exp(median_control))
pValuesSv<-f.pvalue(edata,modSv,mod0Sv)
qValuesSv<-p.adjust(pValuesSv,method="fdr")
table<-cbind(as.data.frame(t(meta_abunds)),as.data.frame(mode),as.data.frame(pValuesSv),as.data.frame(qValuesSv),as.data.frame(log2FC))
colnames(table)<-c("MZ","RT","Mode","PValue","QValue","Log2FC")
POS_RAS<-table%>%filter(Mode=="ESI+")
NEG_RAS<-table%>%filter(Mode=="ESI-")
write.csv(POS_RAS, file="SVA_DE_RAS_POS.csv")
write.csv(NEG_RAS, file="SVA_DE_RAS_NEG.csv")

