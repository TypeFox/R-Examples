### R code from vignette source 'CpGassoc.Rnw'

###################################################
### code chunk number 1: CpGassoc.Rnw:117-136
###################################################
#Sample output from CpGassoc 
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
library(CpGassoc)
data(samplecpg,samplepheno,package="CpGassoc")
results<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)
results
#Analysis with covariates. There are multiple ways to do this. One can define the
#dataframe prior or do it in the function call or as a function such as ~Cov1+Cov2.
#We will do it in the function call
test<-cpg.assoc(samplecpg,samplepheno$weight,data.frame(samplepheno$Distance,samplepheno$Dose),large.data=FALSE)

#Doing a mixed effects model. This does take more time, so we will do a subset of
#the samplecpg
randtest<-cpg.assoc(samplecpg[1:10,],samplepheno$weight,chip.id=samplepheno$chip,random=TRUE,large.data=FALSE)

#summary function will work on items of class cpg.




###################################################
### code chunk number 2: CpGassoc.Rnw:189-197
###################################################
library(CpGassoc)
data(samplecpg,samplepheno,package="CpGassoc")
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
test1<-cpg.assoc(samplecpg[1:100,],samplepheno$weight,large.data=FALSE)
test2<-cpg.assoc(samplecpg[101:200,],samplepheno$weight,large.data=FALSE)
overall<-cpg.combine(list(test1,test2))
overall


###################################################
### code chunk number 3: CpGassoc.Rnw:281-298
###################################################
##Loading the data
library(CpGassoc)
data(samplecpg,samplepheno,package="CpGassoc")
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
#Performing a permutation 10 times
Testperm<-cpg.perm(samplecpg,samplepheno$weight,data.frame(samplepheno$Dose,samplepheno$Distance),
                seed=2314,nperm=10,large.data=FALSE)
Testperm
#All the contents of CpGassoc are included in the output from Testperm
#Using the output from CpGassoc in the example
test<-cpg.assoc(samplecpg,samplepheno$weight,data.frame(samplepheno$Distance,samplepheno$Dose),large.data=FALSE)
all.equal(Testperm$results,test$results)

#summary function works on objects of class cpg.perm
summary(Testperm)



###################################################
### code chunk number 4: CpGassoc.Rnw:350-356
###################################################
library(CpGassoc)
data(samplecpg,samplepheno,package="CpGassoc")
results<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)

cpg.GC(results)
##If the genomic inflation factor is less than one there is no need for adjustment


###################################################
### code chunk number 5: CpGassoc.Rnw:417-418
###################################################
##See the examples in the CpGassoc tutorial.


###################################################
### code chunk number 6: CpGassoc.Rnw:510-512
###################################################
##See the examples listed in cpg.assoc for ways in which to use cpg.work.
##Just change the cpg.assoc to cpg.work.


###################################################
### code chunk number 7: CpGassoc.Rnw:568-581
###################################################
library(CpGassoc)
data(samplecpg,samplepheno,package="CpGassoc")
#Example where there are covariates:
covar<-data.frame(samplepheno$weight,samplepheno$Distance)
test<-design(covar,samplepheno$SBP,samplepheno$chip,FALSE)
dim(test$full)
dim(test$reduced)
test$reduced[1:5,1:5]
test$full[1:5,1:5]
#When no covariates or chip.id:
test2<-design(NULL,samplepheno$SBP,NULL,FALSE)
dim(test2$full)
dim(test2$reduced)


###################################################
### code chunk number 8: manhattan
###################################################
#Doing a Manhattan plot. First load the data:

#Doing a Manhattan plot. First load the data:
library(CpGassoc)
data(samplecpg,samplepheno,annotation,package="CpGassoc")
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
examplemanhat<-cpg.assoc(samplecpg,samplepheno$Disease,large.data=FALSE)

manhattan(examplemanhat,annotation$TargetID,annotation$CHR,annotation$MAPINFO)



###################################################
### code chunk number 9: QQ_plot
###################################################
##Using the results from the example given in cpg.assoc.
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
##QQ Plot:
library(CpGassoc)
data(samplecpg,samplepheno,package="CpGassoc")
test<-cpg.assoc(samplecpg,samplepheno$weight,data.frame(samplepheno$Distance,samplepheno$Dose),large.data=FALSE)
plot(test)
##t-statistic plot:
plot(test,tplot=TRUE)

##Now an example of sort
head(sort(test)$results)

##Summary
summary(test)


###################################################
### code chunk number 10: Perm_qqplot
###################################################
library(CpGassoc)
data(samplecpg,samplepheno,package="CpGassoc")
##We will do the analysis on a subset to save time
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
#The qq plot:
Testperm<-cpg.perm(samplecpg,samplepheno$weight,data.frame(samplepheno$Dose,samplepheno$Distance),
                seed=2314,nperm=10,large.data=FALSE)
plot(Testperm)
#The t-statistic plot from cpg.perm has confidence intervals since we were allowed to perform permutations on the T-values.
plot(Testperm,tplot=TRUE)
#If there was 100 or more permutations, there would be emperical confidence intervals.

###Now for Sort
head(sort(Testperm)$results)
head(Testperm$results)


###################################################
### code chunk number 11: Scatterplot
###################################################
#Load the data:
data(samplecpg,samplepheno,package="CpGassoc")
library(CpGassoc)
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
test<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)
##Using rank, will plot the top three sites in order of significance:
scatterplot(test,c(1:3))
##Using name, specify three sites:
scatterplot(test,cpg.name=c("CpG1182","CpG1000","CpG42"))

##Plotting something that is categorical in nature:
test2<-cpg.assoc(samplecpg,factor(samplepheno$Disease),large.data=FALSE)
scatterplot(test2,c(2))


