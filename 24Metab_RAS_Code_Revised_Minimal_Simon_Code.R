#
###############################################################
######## 24 metabolites panel                                               ###############
###############################################################
#
library(pROC)
library(hlr)
load("metab24.RData")

#gid1 is the group id for the patients in the discovery set, 1: Control, 3: Converter-pre
idx <- c((1:length(gid1))[gid1==1],(1:length(gid1))[gid1==3])

#metabolite level matrix for the patients in the discovery set, row names of the matrix are the metabolite names
matx1 <- matx1[,idx]

dptmp <- gid1[idx]
dptmp1 <- ifelse(dptmp==1,0,1)

#gid2 is the group id for the patients in the validation set, 1: Control, 3: Converter-pre, row names of the matrix are the metabolite names
idx <- c((1:length(gid2))[gid2==1],(1:length(gid2))[gid2==3])

#metabolite level matrix for the patients in the validation set
matx2 <- matx2[,idx]

dptmp <- gid2[idx]
dptmp2 <- ifelse(dptmp==1,0,1)

set.seed(42)
out <- WEMEL(t(matx1),t(matx1),dptmp1,delta=0.001)
phat1 <- coef(out[[2]])[1]+t(matx1)%*%matrix(coef(out[[2]])[-1],nrow(matx1),1)
phat_new = exp(phat1)/(1+exp(phat1))
dy <- factor(dptmp1)

pdf("roc24.pdf",paper="letter",width=8,height=8,bg="transparent")


rocobj.nc.conv.d<- plot.roc(dy, as.vector(phat_new),ylab="True positive rate", xlab="False positive rate", main="ROC curve for the 24 metabolites for the discovery set\n NC vs CONVERTERpre", percent=TRUE,  ci=TRUE,print.auc=TRUE,legacy.axes=TRUE)

ciobj <- ci.se(rocobj.nc.conv.d, specificities=seq(0, 100, 5))
plot(ciobj, type="shape", col="#1c61b6AA")

phat1 <- coef(out[[2]])[1]+t(matx2)%*%matrix(coef(out[[2]])[-1],nrow(matx1),1)
phat_new = exp(phat1)/(1+exp(phat1))
dy <- factor(dptmp2)

rocobj.nc.conv.v<- plot.roc(dy, as.vector(phat_new),ylab="True positive rate", xlab="False positive rate", main="ROC curve for the 24 metabolites for the validation set\n NC vs CONVERTERpre", percent=0.16,  ci=TRUE,print.auc=TRUE,legacy.axes=TRUE)

ciobj <- ci.se(rocobj.nc.conv.v, specificities=seq(0, 100, 5))
plot(ciobj, type="shape", col="#1c61b6AA")

dev.off()

#######SHOULD EVALUATE TRUE
all.equal(rownames(matx1),rownames(matx2))
rownames<-c("Intercept",rownames(matx1))
together<-as.data.frame(cbind(rownames,out$WEMEL))
colnames(together)<-c("Term","Beta Weight")
write.csv(together,file="HLR_Seed42_BetaWeights_RAS_Fiandaca.csv")

save(out,file="24Model_GLM.RData")










