### R code from vignette source 'gskat.Rnw'

###################################################
### code chunk number 1: load.data
###################################################
library(gskat)
data(gdata)
names(gdata)

attach(gdata)
#head(ID);head(y);head(X);head(Z)


###################################################
### code chunk number 2: bintest
###################################################
gskat_score(gdata) #using the kinship working correlation matrix

gskat_score(gdata,F1=FALSE) #using identity working correlation


###################################################
### code chunk number 3: pertu1
###################################################
gskat_score_pert(gdata)$pval 


###################################################
### code chunk number 4: pertu2
###################################################
# gskat_score_pert(gdata,np=100000) #100000 resampling trials



###################################################
### code chunk number 5: pertu3
###################################################
# gskat_score_pert(gdata,pw="Norm") #using normal distribtion instead of Rademacher



###################################################
### code chunk number 6: plink
###################################################
# fileName=foo #plink PED file name
# system(paste("plink --noweb --file ", fileName, " --recodeA --out TEMP/", 
#	fileName,sep="") #plink additive coding
# RAW=read.table(paste("TEMP/",fileName,".raw",sep=""),as.is=T,
#	header=T) #read into R the plink RAW file
# RAW=RAW[order(RAW$FID),] #sort according to Family ID
# ID=RAW[,("FID","IID","PAT","MAT")]
# y=RAW[,"PHENOTYPE"]
# Z=as.matrix(RAW[,-1:-6])
# #X prepared by users
# mydata<-list(ID,y,X,Z)
# gskat_score(mydata) #done


###################################################
### code chunk number 7: sequence1
###################################################
# attach(gdata) #same format as before
# gskat_seq(y,XC=X,Z,ID)
# gskat_seq(y,XC=X,Z,ID,resampling=FALSE) # get asymptoptic p-value only
# gskat_seq(y,XC=X,Z,ID,pw="Norm") #using normal r.v. in pertubation
# gskat_seq(y,XC=X,Z,ID,SNP.weigts)#using customized SNV weighting scores



###################################################
### code chunk number 8: sequence2
###################################################
# attach(gdata) 
# gskat_seq_cont(y,XC=X,Z,ID)
# gskat_seq_cont(y,XC=X,Z,ID, w_a=1,w_b=1) #common variants with no/equal weights


