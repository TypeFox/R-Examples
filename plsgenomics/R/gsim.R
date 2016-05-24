### gsim.R  (2006-01)
###
###                 GSIM for binary data
###
### Copyright 2006-01 Sophie Lambert-Lacroix and Julie Peyre
###
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

gsim <- function(Xtrain, Ytrain, Xtest=NULL, Lambda, hA, hB=NULL, NbIterMax=50)
{
##  INPUT variables
####################
##
##  Xtrain   : matrix ntrain x p
##      train data matrix
##  Ytrain   : vector of length ntrain
##      response variable {1,2}-valued vector
##  Xtest   : matrix ntest x p
##      test data matrix
##  Lambda : real
##      value for the regularization parameter Lambda
##  hA : real
##      indicates h parameter value
##  hB : real
##      indicates h parameter value for second step
##      (if not valid, hB estimated by plug-in)
##  NbIterMax : positive integer
##      max number of iteration in the WIRRLS part
##
##
##  OUTPUT variables
#####################
##  Structure with fields
##      Ytest : vector of length ntest
##      predicted classes for the test data set
##      beta : vector of length p
##      estimated direction for the projection
##      DeletedCol : vector
##           if some covariables have nul variance, DeletedCol gives the 
##               corresponding column number. Otherwise DeletedCol = NULL
##      Cvg : logical
##      1 if convergence in WIRRLS 0 otherwise
##  hB : value of hB used in second step
##
#################################
## TEST on the INPUT variables
#################################
## on X train
if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE))
   stop("Message from gsim.R: Xtrain is not of valid type")

if (dim(Xtrain)[2]==1)
  stop("Message from gsim.R: p=1 is not valid")

ntrain <- dim(Xtrain)[1]
p <- dim(Xtrain)[2]

## on Xtest
ntest <- 0
if (is.null(Xtest)==FALSE)
{
   if (is.vector(Xtest)==TRUE)
      Xtest <- matrix(Xtest,nrow=1)

   if ((is.matrix(Xtest)==FALSE)||(is.numeric(Xtest)==FALSE))
      stop("Message from gsim.R: Xtest is not of valid type")

   if (dim(Xtrain)[2]!=dim(Xtest)[2])
      stop("Message from gsim.R: columns of Xtest and columns of Xtrain must be equal")
  
   ntest <- dim(Xtest)[1] 
}

## on Ytrain

if ((is.vector(Ytrain)==FALSE)||(is.numeric(Ytrain)==FALSE))
   stop("Message from gsim.R: Ytrain is not of valid type")

if (length(Ytrain)!=ntrain)
   stop("Message from gsim.R: the length of Ytrain is not equal to the Xtrain row number")

Ytrain <- Ytrain-1
k0 <- length(which(Ytrain==0))
k1 <- length(which(Ytrain==1))
if (k0+k1 != ntrain)
   stop("Message from gsim.R: Ytrain must be a binary vector (with 0 and 1)")

if ((k0 ==0)||(k1==0))
   stop("Message from gsim.R: there are empty classes")

## on hyper parameters 
# Lambda
if (is.vector(Lambda)==FALSE)
   stop("Message from gsim.R: Lambda is not of valid type")
 
if (length(Lambda)!=1)
   stop("Message from gsim.R: only one value can be specified for Lambda")
   
if ((is.numeric(Lambda)==FALSE)||(Lambda<0))
 stop("Message from gsim.R: Lambda is not of valid type")

# hA
if (is.vector(hA)==FALSE)
   stop("Message from gsim.R: hA is not of valid type")

if (length(hA)!=1)
   stop("Message from gsim.R: only one value can be specified for hA")

if ((is.numeric(hA)==FALSE)||(hA<=0))
 stop("Message from gsim.R: hA is not of valid type")

# hB
if (is.null(hB)==FALSE)
{
   if (is.vector(hB)==FALSE)
      stop("Message from gsim.R: hB is not of valid type1")
   
   if (length(hB)!=1)
      stop("Message from gsim.R: only one value can be specified for hB")

   if ((is.numeric(hB)==FALSE)||(hB<=0))
      stop("Message from gsim.R: hB is not of valid type")
}

# NbIterMax
if (is.vector(NbIterMax)==FALSE)
   stop("Message from gsim.R: NbIterMax is not of valid type")

if ((length(NbIterMax)!=1)||(is.numeric(NbIterMax)==FALSE)||(round(NbIterMax)-NbIterMax!=0)||(NbIterMax<1))
   stop("Message from gsim.R: NbIterMax is not of valid type")

#Some initializations 
r <- min(p,ntrain)
DeletedCol <- NULL
Cvg <- 1
#############################################################
##  STEP 1:
##      1/ Center and Standardize the data matrix
##      2/ Move in the reduced space
##      3/ Form the response variable and the design matrix
#############################################################
#   Center and standardize

MeanXtrain <- apply(Xtrain,2,mean)
sXtrain <- sweep(Xtrain,2,MeanXtrain,FUN="-")
Sigma2train <- apply(sXtrain**2,2,mean)
# Sigma2train <- apply(sXtrain,2,var)*(ntrain-1)/ntrain
if (sum(Sigma2train==0)!=0)
{
   warning("There are covariables with nul variance")
   DeletedCol <- which(Sigma2train==0)
   sXtrain <- sXtrain[,-DeletedCol]
   #names <- names[-DeletedCol]
   
   if (ntest != 0)
      Xtest <- Xtest[,-DeletedCol]
    
   MeanXtrain <-MeanXtrain[-DeletedCol]
   Sigma2train <-Sigma2train[-DeletedCol]
   p <- dim(sXtrain)[2]
   r <- min(p,ntrain)
}

sXtrain <- sweep(sXtrain,2,sqrt(Sigma2train),FUN="/")

# Compute the svd when necessary
if (p>ntrain)
{
   svd.sXtrain <- svd(t(sXtrain))
   D <- diag(c(svd.sXtrain$d))
   r <- length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
   V <- svd.sXtrain$u[,1:r]
   D <- diag(c(svd.sXtrain$d[1:r]))
   U <- svd.sXtrain$v[,1:r]
   sXtrain <- U%*%D
   rm(D)
   rm(U)
   rm(svd.sXtrain)
}

if (ntest != 0)
{
   sXtest <- sweep(Xtest,2,MeanXtrain,FUN="-")
   sXtest <- sweep(sXtest,2,sqrt(Sigma2train),FUN="/")
   
   if (p>ntrain)
   {
      sXtest <- as.matrix(sXtest)%*%V
   }

}

rm(Xtrain)
rm(Xtest)

#   Form the design matrix Ztrain
Ztrain <- cbind(rep(1,ntrain),sXtrain)

#########################################################################
##  STEP 2 :
##      Loop over the nLearn problems
#########################################################################
BETA <- rep(0,length=r)

for (j in 1:ntrain)
{
    #  compute Kernel
    Wk <- exp(-(sweep(sXtrain,2,sXtrain[j,],FUN="-"))**2/(2*hA**2))
    Wk <- apply(Wk,1,prod)
    Wk <- Wk/sum(Wk)
     
    fit <- wirrls(Y=Ytrain,Z=Ztrain,Lambda=Lambda,NbrIterMax=NbIterMax,WKernel=diag(c(Wk)))
    
    if (fit$Cvg == 0)
    {
       warning("Message from gsim.R: the algorithm did not converge in step A")
       Cvg <- 0
    }
    cte <- fit$Coefficients[1,]
    gamma <- fit$Coefficients[-1,]
   
    #constante modification
    cte <- cte+sXtrain[j,]%*%gamma
    expcte <- exp(cte)
    
    gamma <- gamma*expcte/(1+expcte)-(gamma*expcte)*expcte/(1+expcte)^2
    BETA <- BETA + gamma/ntrain

} # end for

rm(gamma)

#########################################################################
##  STEP 3 :
##     hB Estimation step
#########################################################################
# estimating the hB value using plugin
newsXtrain <- sXtrain%*%BETA
if (is.null(hB)==TRUE)
   hB <- (hplugin(newsXtrain,Ytrain,tau=0.1,intgsize=100,NbrIterMax=NbIterMax))$h
   
#########################################################################
##  STEP 4 :
##      Classification step
#########################################################################
# only if Xtest != NULL

if (is.null(Xtest)==FALSE)
{   
   newsXtest <- sXtest%*%BETA 
   ajEst <- rep(1,length=ntest)
   cvg <- 1

   for (kk in 1:ntest)
   {
       Xikdes <- newsXtrain-newsXtest[kk,]
   
       Wk <- exp((-Xikdes**2)/(2*hB^2))
       Xikdes <-cbind(rep(1,length=ntrain),Xikdes)
    
       #new fit over the FDR space
       fit <- wirrls(Ytrain,Xikdes,Lambda=0,NbrIterMax=NbIterMax,WKernel=diag(c(Wk))) 
       
       if (fit$Cvg==0)
       {
          warning("Message from gsim.R: the algorithm did not converge in the step B (i.e. after projection)")
          Cvg <- 0
       }
       ajEst[kk] = fit$Coefficients[1]
    }
}

rm(Xikdes)

#########################################################################
##  STEP 5 :
##      Prediction step
#########################################################################

piest <- 1/(1+exp(-ajEst))
hatY <- as.numeric(piest> 0.5)

#########################################################################
##  STEP 6 :
## CONCLUDE
#########################################################################

##Compute the estimated directions in the initial space

if (p>ntrain)
{beta <- diag(c(1/sqrt(Sigma2train)))%*%V%*%BETA}
if (p<=ntrain)
{beta <- diag(c(1/sqrt(Sigma2train)))%*%BETA}

List <- list(Ytest=(hatY+1),beta=beta,hB=hB,DeletedCol=DeletedCol,Cvg=Cvg)
return(List)

}
