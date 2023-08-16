library(impute)
library(limma)
library(sva)
library(dplyr)


load("RAS_By_Group.RData",verbose=TRUE)
ras_pos<-rbind(converter_POS,control_POS,mciad_POS,super_POS)
ras_neg<-rbind(converter_NEG,control_NEG,mciad_NEG,super_NEG)
all.equal(ras_pos$SampleID,ras_neg$SampleID)

####Read in abundance data and subset to final table of features 
ras_pos<-rbind(converter_POS,control_POS,mciad_POS,super_POS)[,-c(1:2)]
ras_pos_mz<-gsub("_.+$","",colnames(ras_pos))
ras_pos_rt<-gsub("^.+_","",colnames(ras_pos))
ras_neg<-rbind(converter_NEG,control_NEG,mciad_NEG,super_NEG)[,-c(1:2)]
ras_neg_mz<-gsub("_.+$","",colnames(ras_neg))
ras_neg_rt<-gsub("^.+_","",colnames(ras_neg))


meta_pheno<-rbind(converter_POS,control_POS,mciad_POS,super_POS)[,c(1:2)]
meta_abunds<-cbind(rbind(ras_pos_mz,ras_pos_rt),rbind(ras_neg_mz,ras_neg_rt))
abunds<-cbind(ras_pos,ras_neg)
colnames(abunds)<-seq(1,dim(abunds)[2])
mode<-as.factor(c(rep("positive",dim(ras_pos)[2]),rep("negative",dim(ras_neg)[2])))
abunds[abunds==0]<-NA
index<-caret::nearZeroVar(abunds)
meta_abunds<-meta_abunds[,-index]
Mode<-mode[-index]
abunds<-abunds[,-index]

####Impute missing data with KNN imputation
library(impute)
T_abunds<-t(abunds)
edata<-impute.knn(T_abunds,k = 10, rowmax = 1.0, colmax = 1.0, maxp = 1500, rng.seed=362436069)$data
edata<-log2(apply(edata,2,as.numeric))
rownames(edata)<-rownames(T_abunds)

feature_meta<-cbind(as.data.frame(colnames(abunds)),as.data.frame(Mode),as.data.frame(t(meta_abunds)))
feature_meta_use<-cbind(rownames(T_abunds),feature_meta)
colnames(feature_meta_use)[1]<-"Feature"

####Conduct SVA over samples 
pheno_table<-meta_pheno[,1:2]
colnames(pheno_table)<-c("StudyID","DiseaseState")
mod<-model.matrix(~as.factor(DiseaseState),data=pheno_table)
mod0<-model.matrix(~1,data=pheno_table)
n.sv_BE<-num.sv(edata,mod,seed=122)
n.sv_LEEK<-num.sv(edata,mod,seed=122,method="leek")
svobj_LEEK<-sva(edata,mod,mod0,n.sv=n.sv_LEEK)$sv
svobj_BE<-sva(edata,mod,mod0,n.sv=n.sv_BE)$sv
total_LEEK<-cbind(as.data.frame(pheno_table),as.data.frame(svobj_LEEK))
total_BE<-cbind(as.data.frame(pheno_table),as.data.frame(svobj_BE))
colnames(total_LEEK)[1:2]<-c("SampleID","Main")
colnames(total_BE)[1:2]<-c("SampleID","Main")
design_table_BE<-cbind(model.matrix(~0+as.factor(Main),data=total_BE),as.data.frame(svobj_BE))
colnames(design_table_BE)[1:4]<-c("Control","Converter","MCIAD","Super")

design_table_LEEK<-cbind(model.matrix(~0+as.factor(Main),data=total_LEEK),as.data.frame(svobj_LEEK))
colnames(design_table_LEEK)[1:4]<-c("Control","Converter","MCIAD","Super")



####Fit linear model for pairwise contrasts 
arrayw<-arrayWeights(edata, design=design_table_LEEK)  ########design_table_LEEK
fit1<-lmFit(edata,design_table_LEEK,weights=arrayw)
cm1 <- makeContrasts(`Converter-Control`= Converter-Control,
					 levels=design_table_LEEK)
cm2 <- makeContrasts(`Super-Control`= Super-Control,
					 levels=design_table_LEEK)
cm3 <- makeContrasts(`Super-MCIAD`= Super-MCIAD,
					 levels=design_table_LEEK)
cm4 <- makeContrasts(`Super-Converter`= Super-Converter,
					 levels=design_table_LEEK)
cm5 <- makeContrasts(`MCIAD-Control`= MCIAD-Control,
					 levels=design_table_LEEK)
cm6 <- makeContrasts(`MCIAD-Converter`= MCIAD-Converter,
					 levels=design_table_LEEK)
cm6Fix<-  makeContrasts(`Converter-MCIAD`= Converter-MCIAD,
					 levels=design_table_LEEK)


####Fit linear model for pairwise contrasts 
arrayw<-arrayWeights(edata, design=design_table_BE)  ########design_table_BE
fit2<-lmFit(edata,design_table_BE,weights=arrayw)
cm7 <- makeContrasts(`Converter-Control`= Converter-Control,
					 levels=design_table_BE)
cm8 <- makeContrasts(`Super-Control`= Super-Control,
					 levels=design_table_BE)
cm9 <- makeContrasts(`Super-MCIAD`= Super-MCIAD,
					 levels=design_table_BE)
cm10 <- makeContrasts(`Super-Converter`= Super-Converter,
					 levels=design_table_BE)
cm11 <- makeContrasts(`MCIAD-Control`= MCIAD-Control,
					 levels=design_table_BE)
cm12 <- makeContrasts(`MCIAD-Converter`= MCIAD-Converter,
					 levels=design_table_BE)
cm12Fix <- makeContrasts(`Converter-MCIAD`= Converter-MCIAD,
					 levels=design_table_BE)


#####################################PIUmet
#############
#############LEEK
arrayw<-arrayWeights(edata, design=design_table_LEEK)  ########design_table_LEEK
fit1<-lmFit(edata,design_table_LEEK,weights=arrayw)
fit1_F <- contrasts.fit(fit1, cm5)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
LEEK<-joinT%>%select("MZ","Mode","P.Value")






######mutate(`Prize`= -log10(p.adjust(P.Value,method="fdr")))
#############
#############BE
arrayw<-arrayWeights(edata, design=design_table_BE)  ########design_table_BE
fit2<-lmFit(edata,design_table_BE,weights=arrayw)
fit2_F <- contrasts.fit(fit2, cm11)
fit2_F <- eBayes(fit2_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit2_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
BE<-joinT%>%select("MZ","Mode","P.Value")



write.table(rbind(LEEK,BE)%>%mutate(`Prize`= -log10(P.Value))%>%select("MZ","Mode","Prize"),
			sep="\t",file="RAS_PIUmet_5_2022_MCIAD_Control.txt",row.names=FALSE)



#############################Chapter 3
fit1_F <- contrasts.fit(fit1, cm6Fix)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Converter.MCIAD"),
			sep="\t",file="RAS_POS_9_2_Conv_MCIAD_LEEK.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Converter.MCIAD"),
			sep="\t",file="RAS_NEG_9_2_Conv_MCIAD_LEEK.txt",row.names=FALSE)

#############
fit1_F <- contrasts.fit(fit2, cm12Fix)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Converter.MCIAD"),
			sep="\t",file="RAS_POS_9_2_Conv_MCIAD_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Converter.MCIAD"),
			sep="\t",file="RAS_NEG_9_2_Conv_MCIAD_BE.txt",row.names=FALSE)

cd /Volumes/TJGross_Remote/318557/Active_Projects_4_9_2021/ROCAS
mummichog -f RAS_POS_9_2_Conv_MCIAD_LEEK.txt  -o RAS_POS_9_2_Conv_MCIAD_LEEK -m positive 
mummichog -f RAS_NEG_9_2_Conv_MCIAD_LEEK.txt  -o RAS_NEG_9_2_Conv_MCIAD_LEEK -m negative 

mummichog -f RAS_POS_9_2_Conv_MCIAD_BE.txt -o RAS_POS_9_2_Conv_MCIAD_BE -m positive 
mummichog -f RAS_NEG_9_2_Conv_MCIAD_BE.txt -o RAS_NEG_9_2_Conv_MCIAD_BE -m negative












#############
#############
fit1_F <- contrasts.fit(fit1, cm1)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Converter.Control"),
			sep="\t",file="RAS_POS_5_26_Conv_Con_LEEK.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Converter.Control"),
			sep="\t",file="RAS_NEG_5_26_Conv_Con_LEEK.txt",row.names=FALSE)

#############
fit1_F <- contrasts.fit(fit1, cm2)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Super.Control"),
			sep="\t",file="RAS_POS_5_26_Super_Con_LEEK.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Super.Control"),
			sep="\t",file="RAS_NEG_5_26_Super_Con_LEEK.txt",row.names=FALSE)

fit1_F <- contrasts.fit(fit1, cm3)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Super.MCIAD"),
			sep="\t",file="RAS_POS_5_26_Super_MCIAD_LEEK.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Super.MCIAD"),
			sep="\t",file="RAS_NEG_5_26_Super_MCIAD_LEEK.txt",row.names=FALSE)

fit1_F <- contrasts.fit(fit1, cm4)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Super.Converter"),
			sep="\t",file="RAS_POS_5_26_Super_Conv_LEEK.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Super.Converter"),
			sep="\t",file="RAS_NEG_5_26_Super_Conv_LEEK.txt",row.names=FALSE)

fit1_F <- contrasts.fit(fit1, cm5)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","MCIAD.Control"),
			sep="\t",file="RAS_POS_5_26_MCIAD_Con_LEEK.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","MCIAD.Control"),
			sep="\t",file="RAS_NEG_5_26_MCIAD_Con_LEEK.txt",row.names=FALSE)

fit1_F <- contrasts.fit(fit1, cm6)
fit1_F <- eBayes(fit1_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit1_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","MCIAD.Converter"),
			sep="\t",file="RAS_POS_5_26_MCIAD_Conv_LEEK.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","MCIAD.Converter"),
			sep="\t",file="RAS_NEG_5_26_MCIAD_Conv_LEEK.txt",row.names=FALSE)

########LEEK
cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/ROCAS 
mummichog -f RAS_POS_5_26_Conv_Con_LEEK.txt  -o RAS_POS_5_26_Conv_Con_LEEK -m positive 
mummichog -f RAS_NEG_5_26_Conv_Con_LEEK.txt  -o RAS_NEG_5_26_Conv_Con_LEEK -m negative 

mummichog -f RAS_POS_5_26_Super_Con_LEEK.txt -o RAS_POS_5_26_Super_Con_LEEK -m positive 
mummichog -f RAS_NEG_5_26_Super_Con_LEEK.txt -o RAS_NEG_5_26_Super_Con_LEEK -m negative 

mummichog -f RAS_POS_5_26_Super_MCIAD_LEEK.txt -o RAS_POS_5_26_Super_MCIAD_LEEK -m positive 
mummichog -f RAS_NEG_5_26_Super_MCIAD_LEEK.txt -o RAS_NEG_5_26_Super_MCIAD_LEEK -m negative 

mummichog -f RAS_POS_5_26_Super_Conv_LEEK.txt -o RAS_POS_5_26_Super_Conv_LEEK -m positive 
mummichog -f RAS_NEG_5_26_Super_Conv_LEEK.txt -o RAS_NEG_5_26_Super_Conv_LEEK -m negative 

mummichog -f RAS_POS_5_26_MCIAD_Con_LEEK.txt -o RAS_POS_5_26_MCIAD_Con_LEEK -m positive 
mummichog -f RAS_NEG_5_26_MCIAD_Con_LEEK.txt -o RAS_NEG_5_26_MCIAD_Con_LEEK -m negative 

mummichog -f RAS_POS_5_26_MCIAD_Conv_LEEK.txt -o RAS_POS_5_26_MCIAD_Conv_LEEK -m positive 
mummichog -f RAS_NEG_5_26_MCIAD_Conv_LEEK.txt -o RAS_NEG_5_26_MCIAD_Conv_LEEK -m negative 

####################################
####################################
####################################

fit2_F <- contrasts.fit(fit2, cm7)
fit2_F <- eBayes(fit2_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit2_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Converter.Control"),
			sep="\t",file="RAS_POS_5_26_Conv_Con_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Converter.Control"),
			sep="\t",file="RAS_NEG_5_26_Conv_Con_BE.txt",row.names=FALSE)

#############
fit2_F <- contrasts.fit(fit2, cm8)
fit2_F <- eBayes(fit2_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit2_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Super.Control"),
			sep="\t",file="RAS_POS_5_26_Super_Con_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Super.Control"),
			sep="\t",file="RAS_NEG_5_26_Super_Con_BE.txt",row.names=FALSE)

fit2_F <- contrasts.fit(fit2, cm9)
fit2_F <- eBayes(fit2_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit2_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Super.MCIAD"),
			sep="\t",file="RAS_POS_5_26_Super_MCIAD_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Super.MCIAD"),
			sep="\t",file="RAS_NEG_5_26_Super_MCIAD_BE.txt",row.names=FALSE)

fit2_F <- contrasts.fit(fit2, cm10)
fit2_F <- eBayes(fit2_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit2_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","Super.Converter"),
			sep="\t",file="RAS_POS_5_26_Super_Conv_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","Super.Converter"),
			sep="\t",file="RAS_NEG_5_26_Super_Conv_BE.txt",row.names=FALSE)

fit2_F <- contrasts.fit(fit2, cm11)
fit2_F <- eBayes(fit2_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit2_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","MCIAD.Control"),
			sep="\t",file="RAS_POS_5_26_MCIAD_Con_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","MCIAD.Control"),
			sep="\t",file="RAS_NEG_5_26_MCIAD_Con_BE.txt",row.names=FALSE)

fit2_F <- contrasts.fit(fit2, cm12)
fit2_F <- eBayes(fit2_F,trend=TRUE)

####Output final DE table and write to file
T<-topTableF(fit2_F,adjust="BH",number=100000)
T<-cbind(rownames(T),T)
colnames(T)[1]<-"Feature"
joinT<-dplyr::inner_join(feature_meta_use,T)
colnames(joinT)[c(4,5)]<-c("MZ","RT")
write.table(joinT%>%filter(Mode=="ESI+")%>%select("MZ","RT","P.Value","MCIAD.Converter"),
			sep="\t",file="RAS_POS_5_26_MCIAD_Conv_BE.txt",row.names=FALSE)
write.table(joinT%>%filter(Mode=="ESI-")%>%select("MZ","RT","P.Value","MCIAD.Converter"),
			sep="\t",file="RAS_NEG_5_26_MCIAD_Conv_BE.txt",row.names=FALSE)

########BE
cd /Volumes/NO\ NAME/Active_Projects_4_9_2021/ROCAS 
mummichog -f RAS_POS_5_26_Conv_Con_BE.txt  -o RAS_POS_5_26_Conv_Con_BE -m positive 
mummichog -f RAS_NEG_5_26_Conv_Con_BE.txt  -o RAS_NEG_5_26_Conv_Con_BE -m negative 

mummichog -f RAS_POS_5_26_Super_Con_BE.txt -o RAS_POS_5_26_Super_Con_BE -m positive 
mummichog -f RAS_NEG_5_26_Super_Con_BE.txt -o RAS_NEG_5_26_Super_Con_BE -m negative 

mummichog -f RAS_POS_5_26_Super_MCIAD_BE.txt -o RAS_POS_5_26_Super_MCIAD_BE -m positive 
mummichog -f RAS_NEG_5_26_Super_MCIAD_BE.txt -o RAS_NEG_5_26_Super_MCIAD_BE -m negative 

mummichog -f RAS_POS_5_26_Super_Conv_BE.txt -o RAS_POS_5_26_Super_Conv_BE -m positive 
mummichog -f RAS_NEG_5_26_Super_Conv_BE.txt -o RAS_NEG_5_26_Super_Conv_BE -m negative 

mummichog -f RAS_POS_5_26_MCIAD_Con_BE.txt -o RAS_POS_5_26_MCIAD_Con_BE -m positive 
mummichog -f RAS_NEG_5_26_MCIAD_Con_BE.txt -o RAS_NEG_5_26_MCIAD_Con_BE -m negative 

mummichog -f RAS_POS_5_26_MCIAD_Conv_BE.txt -o RAS_POS_5_26_MCIAD_Conv_BE -m positive 
mummichog -f RAS_NEG_5_26_MCIAD_Conv_BE.txt -o RAS_NEG_5_26_MCIAD_Conv_BE -m negative 
