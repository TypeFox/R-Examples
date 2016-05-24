pkgname <- "plsgenomics"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('plsgenomics')

assign(".oldSearch", search(), pos = 'CheckExEnv')

###############################################################################################
cleanEx()
nameEx("Colon")
### * Colon

flush(stderr()); flush(stdout())

### Name: Colon
### Title: Gene expression data from Alon et al. (1999)
### Aliases: Colon
### Keywords: datasets

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load data set
data(Colon)

# how many samples and how many genes ?
dim(Colon$X)

# how many samples of class 1 and 2 respectively ?
sum(Colon$Y==1)
sum(Colon$Y==2)


###############################################################################################
cleanEx()
nameEx("Ecoli")
### * Ecoli

flush(stderr()); flush(stdout())

### Name: Ecoli
### Title: Ecoli gene expression and connectivity data from Kao et al.
###   (2003)
### Aliases: Ecoli
### Keywords: datasets

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load data set
data(Ecoli)

# how many genes and how many transcription factors ?
dim(Ecoli$CONNECdata)




###############################################################################################
cleanEx()
nameEx("SRBCT")
### * SRBCT

flush(stderr()); flush(stdout())

### Name: SRBCT
### Title: Gene expression data from Khan et al. (2001)
### Aliases: SRBCT
### Keywords: datasets

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load data set
data(SRBCT)

# how many samples and how many genes ?
dim(SRBCT$X)

# how many samples of class 1, 2, 3 and 4, respectively ?
sum(SRBCT$Y==1)
sum(SRBCT$Y==2)
sum(SRBCT$Y==3)
sum(SRBCT$Y==4)



###############################################################################################
cleanEx()
nameEx("TFA.estimate")
### * TFA.estimate

flush(stderr()); flush(stdout())

### Name: TFA.estimate
### Title: Prediction of Transcription Factor Activities using PLS
### Aliases: TFA.estimate
### Keywords: regression

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load Ecoli data
data(Ecoli)

# estimate TFAs based on 3 latent components
TFA.estimate(Ecoli$CONNECdata,Ecoli$GEdata,ncomp=3,nruncv=0)

# estimate TFAs and determine the best number of latent components simultaneously
TFA.estimate(Ecoli$CONNECdata,Ecoli$GEdata,ncomp=1:5,nruncv=20)



###############################################################################################
cleanEx()
nameEx("gsim")
### * gsim

flush(stderr()); flush(stdout())

### Name: gsim
### Title: GSIM for binary data
### Aliases: gsim

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),12),sample(which(Colon$Y==1),8))

Xtrain <- Colon$X[IndexLearn,]
Ytrain <- Colon$Y[IndexLearn]
Xtest <- Colon$X[-IndexLearn,]

# preprocess data
resP <- preprocess(Xtrain= Xtrain, Xtest=Xtest,Threshold = c(100,16000),Filtering=c(5,500),log10.scale=TRUE,row.stand=TRUE)

# perform prediction by GSIM
res <- gsim(Xtrain=resP$pXtrain,Ytrain= Ytrain,Xtest=resP$pXtest,Lambda=10,hA=50,hB=NULL)
   
res$Cvg
sum(res$Ytest!=Colon$Y[-IndexLearn])



###############################################################################################
cleanEx()
nameEx("gsim.cv")
### * gsim.cv

flush(stderr()); flush(stdout())

### Name: gsim.cv
### Title: Determination of the ridge regularization parameter and the
###   bandwidth to be used for classification with GSIM for binary data
### Aliases: gsim.cv

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),12),sample(which(Colon$Y==1),8))

Xtrain <- Colon$X[IndexLearn,]
Ytrain <- Colon$Y[IndexLearn]
Xtest <- Colon$X[-IndexLearn,]

# preprocess data
resP <- preprocess(Xtrain= Xtrain, Xtest=Xtest,Threshold = c(100,16000),Filtering=c(5,500),log10.scale=TRUE,row.stand=TRUE)

# Determine optimum h and lambda
hl <- gsim.cv(Xtrain=resP$pXtrain,Ytrain=Ytrain,hARange=c(7,20),LambdaRange=c(0.1,1),hB=NULL)

# perform prediction by GSIM  
res <- gsim(Xtrain=resP$pXtrain,Ytrain=Ytrain,Xtest=resP$pXtest,Lambda=hl$Lambda,hA=hl$hA,hB=NULL)
res$Cvg
sum(res$Ytest!=Colon$Y[-IndexLearn])



###############################################################################################
cleanEx()
nameEx("leukemia")
### * leukemia

flush(stderr()); flush(stdout())

### Name: leukemia
### Title: Gene expression data from Golub et al. (1999)
### Aliases: leukemia
### Keywords: datasets

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load data set
data(leukemia)

# how many samples and how many genes ?
dim(leukemia$X)

# how many samples of class 1 and 2, respectively ?
sum(leukemia$Y==1)
sum(leukemia$Y==2)


###############################################################################################
cleanEx()
nameEx("mgsim")
### * mgsim

flush(stderr()); flush(stdout())

### Name: mgsim
### Title: GSIM for categorical data
### Aliases: mgsim

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load SRBCT data
data(SRBCT)
IndexLearn <- c(sample(which(SRBCT$Y==1),10),sample(which(SRBCT$Y==2),4),sample(which(SRBCT$Y==3),7),sample(which(SRBCT$Y==4),9))

# perform prediction by MGSIM
res <- mgsim(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],Lambda=0.001,h=19,Xtest=SRBCT$X[-IndexLearn,])
res$Cvg
sum(res$Ytest!=SRBCT$Y[-IndexLearn])

# prediction for another sample
Xnew <- SRBCT$X[83,]
# projection of Xnew onto the c estimated direction
Xproj <- Xnew %*% res$beta
# Compute the linear predictor for each classes expect class 1
eta <- diag(cbind(rep(1,3),t(Xproj)) %*% res$Coefficients)
Ypred <- which.max(c(0,eta))
Ypred
SRBCT$Y[83]



###############################################################################################
cleanEx()
nameEx("mgsim.cv")
### * mgsim.cv

flush(stderr()); flush(stdout())

### Name: mgsim.cv
### Title: Determination of the ridge regularization parameter and the
###   bandwidth to be used for classification with GSIM for categorical
###   data
### Aliases: mgsim.cv

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load SRBCT data
data(SRBCT)
IndexLearn <- c(sample(which(SRBCT$Y==1),10),sample(which(SRBCT$Y==2),4),
                sample(which(SRBCT$Y==3),7),sample(which(SRBCT$Y==4),9))

### Determine optimum h and lambda
# /!\ take 30 secondes to run
#hl <- mgsim.cv(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],
#                            LambdaRange=c(0.1),hRange=c(7,20))

### perform prediction by MGSIM
#res <- mgsim(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],Lambda=hl$Lambda,
#             h=hl$h,Xtest=SRBCT$X[-IndexLearn,])
#res$Cvg
#sum(res$Ytest!=SRBCT$Y[-IndexLearn])



###############################################################################################
cleanEx()
nameEx("mrpls")
### * mrpls

flush(stderr()); flush(stdout())

### Name: mrpls
### Title: Ridge Partial Least Square for categorical data
### Aliases: mrpls

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load SRBCT data
data(SRBCT)
IndexLearn <- c(sample(which(SRBCT$Y==1),10),sample(which(SRBCT$Y==2),4),sample(which(SRBCT$Y==3),7),sample(which(SRBCT$Y==4),9))

# perform prediction by MRPLS
res <- mrpls(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],Lambda=0.001,ncomp=2,Xtest=SRBCT$X[-IndexLearn,])
sum(res$Ytest!=SRBCT$Y[-IndexLearn])

# prediction for another sample
Xnew <- SRBCT$X[83,]
# Compute the linear predictor for each classes expect class 1
eta <- diag(t(cbind(c(1,Xnew),c(1,Xnew),c(1,Xnew))) %*% res$Coefficients)
Ypred <- which.max(c(0,eta))
Ypred
SRBCT$Y[83]



###############################################################################################
cleanEx()
nameEx("mrpls.cv")
### * mrpls.cv

flush(stderr()); flush(stdout())

### Name: mrpls.cv
### Title: Determination of the ridge regularization parameter and the
###   number of PLS components to be used for classification with RPLS for
###   categorical data
### Aliases: mrpls.cv

### ** Examples

# load plsgenomics library
# load plsgenomics library
library(plsgenomics)

# load SRBCT data
data(SRBCT)
IndexLearn <- c(sample(which(SRBCT$Y==1),10),sample(which(SRBCT$Y==2),4),sample(which(SRBCT$Y==3),7),sample(which(SRBCT$Y==4),9))

# Determine optimum ncomp and Lambda
nl <- mrpls.cv(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],LambdaRange=c(0.1,1),ncompMax=3)

# perform prediction by MRPLS
res <- mrpls(Ytrain=SRBCT$Y[IndexLearn],Xtrain=SRBCT$X[IndexLearn,],Lambda=nl$Lambda,ncomp=nl$ncomp,Xtest=SRBCT$X[-IndexLearn,])
sum(res$Ytest!=SRBCT$Y[-IndexLearn])



###############################################################################################
cleanEx()
nameEx("pls.lda")
### * pls.lda

flush(stderr()); flush(stdout())

### Name: pls.lda
### Title: Classification with PLS Dimension Reduction and Linear
###   Discriminant Analysis
### Aliases: pls.lda
### Keywords: multivariate

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load leukemia data
data(leukemia)

# Classify observations 1,2,3 (test set) using observations 4 to 38 (training set), with 2 PLS components
pls.lda(Xtrain=leukemia$X[-(1:3),],Ytrain=leukemia$Y[-(1:3)],Xtest=leukemia$X[1:3,],ncomp=2,nruncv=0)

# Classify observations 1,2,3 (test set) using observations 4 to 38 (training set), with the best number of components as determined by cross-validation
pls.lda(Xtrain=leukemia$X[-(1:3),],Ytrain=leukemia$Y[-(1:3)],Xtest=leukemia$X[1:3,],ncomp=1:4,nruncv=20)



###############################################################################################
cleanEx()
nameEx("pls.lda.cv")
### * pls.lda.cv

flush(stderr()); flush(stdout())

### Name: pls.lda.cv
### Title: Determination of the number of latent components to be used for
###   classification with PLS and LDA
### Aliases: pls.lda.cv
### Keywords: multivariate

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load leukemia data
data(leukemia)

# Determine the best number of components to be used for classification using the cross-validation procedure
# choose the best number from 2,3,4
pls.lda.cv(Xtrain=leukemia$X,Ytrain=leukemia$Y,ncomp=2:4,nruncv=20)
# choose the best number from 1,2,3
pls.lda.cv(Xtrain=leukemia$X,Ytrain=leukemia$Y,ncomp=3,nruncv=20)




###############################################################################################
cleanEx()
nameEx("pls.regression")
### * pls.regression

flush(stderr()); flush(stdout())

### Name: pls.regression
### Title: Multivariate Partial Least Squares Regression
### Aliases: pls.regression
### Keywords: multivariate

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load the Ecoli data
data(Ecoli)

# perform pls regression
# with unit latent components
pls.regression(Xtrain=Ecoli$CONNECdata,Ytrain=Ecoli$GEdata,Xtest=Ecoli$CONNECdata,ncomp=1:3,unit.weights=FALSE)

# with unit weight vectors
pls.regression(Xtrain=Ecoli$CONNECdata,Ytrain=Ecoli$GEdata,Xtest=Ecoli$CONNECdata,ncomp=1:3,unit.weights=TRUE)




###############################################################################################
cleanEx()
nameEx("pls.regression.cv")
### * pls.regression.cv

flush(stderr()); flush(stdout())

### Name: pls.regression.cv
### Title: Determination of the number of latent components to be used in
###   PLS regression
### Aliases: pls.regression.cv
### Keywords: multivariate

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load Ecoli data
data(Ecoli)

# determine the best number of components for PLS regression using the cross-validation approach
# choose the best number from 1,2,3,4
pls.regression.cv(Xtrain=Ecoli$CONNECdata,Ytrain=Ecoli$GEdata,ncomp=4,nruncv=20)
# choose the best number from 2,3
pls.regression.cv(Xtrain=Ecoli$CONNECdata,Ytrain=Ecoli$GEdata,ncomp=c(2,3),nruncv=20)



###############################################################################################
cleanEx()
nameEx("preprocess")
### * preprocess

flush(stderr()); flush(stdout())

### Name: preprocess
### Title: preprocess for microarray data
### Aliases: preprocess

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),27),sample(which(Colon$Y==1),14))

Xtrain <- Colon$X[IndexLearn,]
Ytrain <- Colon$Y[IndexLearn]
Xtest <- Colon$X[-IndexLearn,]

# preprocess data
resP <- preprocess(Xtrain= Xtrain, Xtest=Xtest,Threshold = c(100,16000),Filtering=c(5,500),log10.scale=TRUE,row.stand=TRUE)

# how many genes after preprocess ?
dim(resP$pXtrain)[2]




###############################################################################################
cleanEx()
nameEx("rpls")
### * rpls

flush(stderr()); flush(stdout())

### Name: rpls
### Title: Ridge Partial Least Square for binary data
### Aliases: rpls

### ** Examples

# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),12),sample(which(Colon$Y==1),8))

# preprocess data
res <- preprocess(Xtrain= Colon$X[IndexLearn,], Xtest=Colon$X[-IndexLearn,],Threshold = c(100,16000),Filtering=c(5,500),log10.scale=TRUE,row.stand=TRUE)
# the results are given in res$pXtrain and res$pXtest

# perform prediction by RPLS
resrpls <- rpls(Ytrain=Colon$Y[IndexLearn],Xtrain=res$pXtrain,Lambda=0.6,ncomp=1,Xtest=res$pXtest)
resrpls$hatY
sum(resrpls$Ytest!=Colon$Y[-IndexLearn])

# prediction for another sample
Xnew <- res$pXtest[1,]
# Compute the linear predictor for each classes expect class 0
eta <- c(1,Xnew) %*% resrpls$Coefficients
Ypred <- which.max(c(0,eta))
Ypred





###############################################################################################
cleanEx()
nameEx("rpls.cv")
### * rpls.cv

flush(stderr()); flush(stdout())

### Name: rpls.cv
### Title: Determination of the ridge regularization parameter and the
###   number of PLS components to be used for classification with RPLS for
###   binary data
### Aliases: rpls.cv

### ** Examples

# load plsgenomics library
# load plsgenomics library
library(plsgenomics)

# load Colon data
data(Colon)
IndexLearn <- c(sample(which(Colon$Y==2),12),sample(which(Colon$Y==1),8))

# preprocess data
res <- preprocess(Xtrain= Colon$X[IndexLearn,], Xtest=Colon$X[-IndexLearn,],Threshold = c(100,16000),Filtering=c(5,500),log10.scale=TRUE,row.stand=TRUE)
# the results are given in res$pXtrain and res$pXtest

# Determine optimum ncomp and lambda
nl <- rpls.cv(Ytrain=Colon$Y[IndexLearn],Xtrain=res$pXtrain,LambdaRange=c(0.1,1),ncompMax=3)

# perform prediction by RPLS
resrpls <- rpls(Ytrain=Colon$Y[IndexLearn],Xtrain=res$pXtrain,Lambda=nl$Lambda,ncomp=nl$ncomp,Xtest=res$pXtest)
sum(resrpls$Ytest!=Colon$Y[-IndexLearn])



###############################################################################################
cleanEx()
nameEx("variable.selection")
### * variable.selection

flush(stderr()); flush(stdout())

### Name: variable.selection
### Title: Variable selection using the PLS weights
### Aliases: variable.selection
### Keywords: multivariate

### ** Examples

# load plsgenomics library
library(plsgenomics)

# generate X and Y (4 observations and 3 variables)
X<-matrix(c(4,3,3,4,1,0,6,7,3,5,5,9),4,3,byrow=FALSE)
Y<-c(1,1,2,2)

# select the 2 best variables
variable.selection(X,Y,nvar=2)
# order the 3 variables
variable.selection(X,Y)

# load the leukemia data 
data(leukemia)

# select the 50 best variables from the leukemia data
variable.selection(leukemia$X,leukemia$Y,nvar=50)




###############################################################################################
cleanEx()
nameEx("rirls.spls")
### * rirls.spls

flush(stderr()); flush(stdout())

### Name: rirls.spls
### Title: Classification by Ridge Iteratively Reweighted Least Squares 
###        followed by Adaptive Sparse PLS regression for binary response
### Aliases: rirls.spls
### Keywords: multivariate

### ** Examples

### load plsgenomics library
library(plsgenomics)

### generating data
n <- 50
p <- 100
sample1 <- sample.bin(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, 
                      mean.H=0.2, sigma.H=10, sigma.F=5)

X <- sample1$X
Y <- sample1$Y

### splitting between learning and testing set
index.train <- sort(sample(1:n, size=round(0.7*n)))
index.test <- (1:n)[-index.train]

Xtrain <- X[index.train,]
Ytrain <- Y[index.train,]

Xtest <- X[index.test,]
Ytest <- Y[index.test,]

### fitting the model, and predicting new observations
model1 <- rirls.spls(Xtrain=Xtrain, Ytrain=Ytrain, lambda.ridge=2, lambda.l1=0.5, ncomp=2, 
                     Xtest=Xtest, adapt=TRUE, maxIter=100, svd.decompose=TRUE)
str(model1)

### prediction error rate
sum(model1$hatYtest!=Ytest) / length(index.test)





###############################################################################################
cleanEx()
nameEx("rirls.spls.tune")
### * rirls.spls.tune

flush(stderr()); flush(stdout())

### Name: rirls.spls.tune
### Title: Tuning parameters (ncomp, lambda.l1, lambda.ridge) for Ridge Iteratively Reweighted Least Squares 
###        followed by Adaptive Sparse PLS regression for binary response, by K-fold cross-validation
### Aliases: rirls.spls.tune
### Keywords: multivariate

### ** Examples

### load plsgenomics library
library(plsgenomics)

### generating data
n <- 50
p <- 100
sample1 <- sample.bin(n=n, p=p, kstar=20, lstar=2, beta.min=0.25, beta.max=0.75, mean.H=0.2, 
                      sigma.H=10, sigma.F=5)

X <- sample1$X
Y <- sample1$Y

### hyper-parameters values to test
lambda.l1.range <- seq(0.05,0.95,by=0.3) # between 0 and 1
ncomp.range <- 1:2

# log-linear range between 0.01 a,d 1000 for lambda.ridge.range
logspace <- function( d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 
lambda.ridge.range <- signif(logspace(d1 <- -2, d2 <- 3, n=6), digits=3)

### tuning the hyper-parameters
cv1 <- rirls.spls.tune(X=X, Y=Y, lambda.ridge.range=lambda.ridge.range, 
                       lambda.l1.range=lambda.l1.range, ncomp.range=ncomp.range, 
                       adapt=TRUE, maxIter=100, svd.decompose=TRUE, 
                       return.grid=TRUE, ncores=1, nfolds=10)
str(cv1)







###############################################################################################

### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
