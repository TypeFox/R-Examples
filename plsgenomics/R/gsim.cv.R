### gsim.cv.R  (2006-02)
###
###    Determination by Cross-validation of GSIM
###                 hyper-parameters for binary data
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

gsim.cv <- function(Xtrain, Ytrain,LambdaRange,hARange,hB=NULL,NbIterMax=50)
{
##  INPUT variables
####################
##
##  Xtrain   : matrix ntrain x p
##      train data matrix
##  Ytrain   : vector of length ntrain
##      response variable {1,2}-valued vector
##  LambdaRange : vector nLambda
##      possible values for the regularization parameter Lambda
##  hRandge : vector nh
##      possible values for the bandwidth parameter
##  hB : real
##      indicates h parameter value for second step
##      (if not valid, hB estimated by plug-in)
##  NbIterMax : positive integer
##      max number of iteration in the WIRRLS part
##
##    OUTPUT VARIABLES
##########################
## Lambda : optimal regularization parameter Lambda
## hA : optimal bandwidth parameter

#################################
## TEST on the INPUT variables
################################# 
## on X train
if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE))
   stop("Message from gsim.cv.R: Xtrain is not of valid type")

if (dim(Xtrain)[2]==1)
  stop("Message from gsim.cv.R: p=1 is not valid")

ntrain <- dim(Xtrain)[1]
p <- dim(Xtrain)[2]

## on Ytrain

if ((is.vector(Ytrain)==FALSE)||(is.numeric(Ytrain)==FALSE))
   stop("Message from gsim.cv.R: Ytrain is not of valid type")

if (length(Ytrain)!=ntrain)
   stop("Message from gsim.cv.R: the length of Ytrain is not equal to the Xtrain row number")

Ytrain <- Ytrain-1
k0 <- length(which(Ytrain==0))
k1 <- length(which(Ytrain==1))
if (k0+k1 != ntrain)
   stop("Message from gsim.cv.R: Ytrain must be a binary vector (with 0 and 1)")

if ((k0 <=1)||(k1<=1))
   stop("Message from gsim.cv.R: there are not enough samples for each class")

## on hyper parameters range

if ((is.numeric(LambdaRange)==FALSE)||(is.vector(LambdaRange)==FALSE)||(sum(LambdaRange<0)>0))
   stop("Message from gsim.cv.R: LambdaRange is not of valid type")

LambdaRange <- sort(LambdaRange)

if ((is.numeric(hARange)==FALSE)||(is.vector(hARange)==FALSE)||(sum(hARange<=0)>0))
   stop("Message from gsim.cv.R: hARange is not of valid type")

hARange <- sort(hARange)
 
## on hB
if (is.null(hB)==FALSE)
{
   if (is.vector(hB)==FALSE)
      stop("Message from gsim.cv.R: hB is not of valid type1")
   
   if (length(hB)!=1)
      stop("Message from gsim.cv.R: only one value can be specified for hB")

   if ((is.numeric(hB)==FALSE)||(hB<=0))
      stop("Message from gsim.cv.R: hB is not of valid type")
}

## NbIterMax
if ((is.numeric(NbIterMax)==FALSE)||(round(NbIterMax)-NbIterMax!=0)||(NbIterMax<1))
   stop("Message from gsim.cv.R: NbIterMax is not of valid type")

## CV LOOP
############
#Some initializations

ntrainCV <- ntrain - 1
nL <- length(LambdaRange)
nh <- length(hARange)

ResCV <- matrix(0,nrow=nL,ncol=nh)
CVG <- matrix(1,nrow=nL,ncol=nh)

for (ncv in 1:ntrain)
{
   # determine the data matrix (without ncvth sample)
   cvXtrain <- Xtrain[-ncv,]
   cvXtest <- matrix(Xtrain[ncv,],nrow=1)

   p <- dim(cvXtrain)[2]
   r <- min(p,ntrainCV)

   #   Center and standardize

   MeancvXtrain <- apply(cvXtrain,2,mean)
   sXtrain <- sweep(cvXtrain,2,MeancvXtrain,FUN="-")
   Sigma2train <- apply(sXtrain**2,2,mean)

   if (sum(Sigma2train==0)!=0)
   {
      DeletedCol <- which(Sigma2train==0)
      sXtrain <- sXtrain[,-DeletedCol]
      p <- dim(sXtrain)[2]
      r <- min(p,ntrainCV)
     
      MeancvXtrain <-MeancvXtrain[-DeletedCol]
      Sigma2train <-Sigma2train[-DeletedCol]
      
      rm(DeletedCol)
   }

   sXtrain <- sweep(sXtrain,2,sqrt(Sigma2train),FUN="/")   
   
   # move in the reduced space when necessary
   if (p>ntrainCV)
   {
      svd.sXtrain <- svd(t(sXtrain))
      r<-length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
      V <- svd.sXtrain$u[,1:r]
      D <- diag(c(svd.sXtrain$d[1:r]))
      U <- svd.sXtrain$v[,1:r]
      sXtrain <- U%*%D
      rm(D)
      rm(U)
      rm(svd.sXtrain)
   }

   sXtest <- sweep(cvXtest,2,MeancvXtrain,FUN="-")
   sXtest <- cvXtest - MeancvXtrain
   sXtest <- sweep(sXtest,2,sqrt(Sigma2train),FUN="/")

   if (p>ntrainCV)
   {
      sXtest <- matrix(sXtest,nrow=1)
      sXtest <- sXtest%*%V
      rm(V)
   }
   
   rm(cvXtrain)
   
   #   Form the design matrix Ztrain
   Ztrain <- cbind(rep(1,ntrainCV),sXtrain)
   Ztest <- cbind(rep(1,1),sXtest)
   
   for (j in 1:nh)
   {
      # compute the Kernel matrix
      WKernel <- matrix(0,ntrainCV,ntrainCV)
      for (jj in 1:ntrainCV)
      {
         WKernel[,jj] <- apply(exp(-(sweep(sXtrain,2,sXtrain[jj,],FUN="-"))**2/(2*hARange[j]**2)),1,prod)
         WKernel[,jj] <- WKernel[,jj]/sum(WKernel[,jj])
      }
      
      for (i in 1:nL)
      {

         if (CVG[i,j]==1)
     {
            res <- gsim.aux(Ztrain=Ztrain,Ytrain=Ytrain[-ncv],Ztest=Ztest,Lambda=LambdaRange[i],WKernel,NbIterMax=NbIterMax,hB=hB)

        if (res$Convergence==0)
            {
           CVG[i,j] <- 0
               ResCV[i,j] <- ntrain
        }
    
            if (res$Convergence==1)
            {
           ResCV[i,j] <- ResCV[i,j]+abs(res$hatY-Ytrain[ncv])
        }
     } # if
      } # for i
   } # for j

} # for ncv

## CONCLUDE
##############

if (sum(CVG)==0)
{stop("Message from gsim.cv.R : no optimal Lambda and h for the given Range")}

#else
# Determine optimal Lambda and h
aux <- which.min(ResCV)
if (nh==1)
{
   Lambda <- LambdaRange[aux]
   h <- hARange[1]
}

if (nh!=1)
{
   h <- hARange[(aux-1)%/%nL+1]
   Lambda <- LambdaRange[(aux-1)%%nL+1]
}

return(list(Lambda=Lambda,hA=h))

}

