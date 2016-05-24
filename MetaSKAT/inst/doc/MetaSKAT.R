### R code from vignette source 'MetaSKAT.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: data
###################################################

library(MetaSKAT)
data(Example)
names(Example)
attach(Example)


###################################################
### code chunk number 2: MetaSKAT1
###################################################
# continuous trait 
obj<-Meta_Null_Model(y.list, x.list, n.cohort=3, out_type="D")


# rho=0 (SKAT)
MetaSKAT_wZ(Z.list[[1]], obj)$p.value

# rho=1 (weighted burden test)
MetaSKAT_wZ(Z.list[[1]], obj, r.corr=1)$p.value


# SKAT-O
MetaSKAT_wZ(Z.list[[1]], obj, method="optimal")$p.value



###################################################
### code chunk number 3: MetaSKAT11
###################################################
# rho=0 (SKAT)
MetaSKAT_wZ(Z.list[[1]], obj, is.separate = TRUE, combined.weight=FALSE )$p.value

# SKAT-O
MetaSKAT_wZ(Z.list[[1]], obj, method="optimal", is.separate = TRUE, 
combined.weight=FALSE)$p.value



###################################################
### code chunk number 4: MetaSKAT11
###################################################
# rho=0 (SKAT). First two cohorts belong to the same group
MetaSKAT_wZ(Z.list[[1]], obj, is.separate = TRUE
, combined.weight=FALSE, Group_Idx=c(1,1,2))$p.value

# SKAT-O. First two cohorts belong to the same group
MetaSKAT_wZ(Z.list[[1]], obj, method="optimal"
, is.separate = TRUE, combined.weight=FALSE, Group_Idx=c(1,1,2))$p.value



###################################################
### code chunk number 5: MetaSKAT2
###################################################

File.SetID<-"./01.SetID"
File.Bed<-"./01.bed"
File.Bim<-"./01.bim"
File.Fam<-"./01.fam"
File.Mat<-"./01.MSSD"
File.SetInfo<-"./01.MInfo"

		
FAM<-read.table(File.Fam, header=FALSE)
y<-FAM[,6]

#########################################
# Test Main File

# need SKAT package to use SKAT_Null_Model function
library(SKAT)

N.Sample<-length(y)
obj<-SKAT_Null_Model(y~1)
Generate_Meta_Files(obj, File.Bed, File.Bim
, File.SetID, File.Mat, File.SetInfo,N.Sample)
		


###################################################
### code chunk number 6: MetaSKAT3
###################################################

for( IDX_G in 2:3){

	File.SetID<-sprintf("./%02d.SetID",IDX_G)
	File.Bed<-sprintf("./%02d.bed",IDX_G)
	File.Bim<-sprintf("./%02d.bim",IDX_G)
	File.Fam<-sprintf("./%02d.fam",IDX_G)
	File.Mat<-sprintf("./%02d.MSSD",IDX_G)
	File.SetInfo<-sprintf("./%02d.MInfo",IDX_G)
		
	FAM<-read.table(File.Fam, header=FALSE)
	y<-FAM[,6]
	N.Sample<-length(y)
	obj<-SKAT_Null_Model(y~1)
	re1<-Generate_Meta_Files(obj, File.Bed, File.Bim, 
	File.SetID, File.Mat, File.SetInfo, N.Sample)
	
}



###################################################
### code chunk number 7: MetaSKAT4
###################################################

File.Mat.vec<-rep("",3)
File.Info.vec<-rep("",3)

for( IDX_G in 1:3){

	File.Mat<-sprintf("./%02d.MSSD",IDX_G)
	File.Info<-sprintf("./%02d.MInfo",IDX_G)
		
	File.Mat.vec[IDX_G]<-File.Mat
	File.Info.vec[IDX_G]<-File.Info
		
}

# open files

Cohort.Info<-Open_MSSD_File_2Read(File.Mat.vec, File.Info.vec)

# get a p-value of the first set.

MetaSKAT_MSSD_OneSet(Cohort.Info, SetID="1")$p.value

# get p-values of all sets

MetaSKAT_MSSD_ALL(Cohort.Info)



