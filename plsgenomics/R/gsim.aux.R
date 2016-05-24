### gsim.aux.R  (2006-01)
###
###                 A GSIM auxiliary function for Cross Validation
###             for binary data
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

gsim.aux <- function(Ztrain, Ytrain, Ztest, Lambda,WKernel, hB=NULL, NbIterMax=15)
{
##  INPUT variables
####################
##
##  Ztrain   : matrix ntrain x (r+1)
##      train design matrix
##  Ytrain   : vector of length ntrain
##      response variable {0,1}-valued vector
##  Ztest   : matrix 1 x r
##      test design matrix
##  Lambda : real
##      value for the regularization parameter Lambda
##  WKernel : matrix ntrain x ntrain
##      bandwith matrix
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
##     hatY : vector of length ntest
##      predicted classes for the test data set
##      Convergence : logical
##      did the algorithm converge ?

ntrain <- dim(Ztrain)[1]
r <- dim(Ztrain)[2] - 1

if (is.vector(Ztest)==TRUE)
   Ztest <- matrix(Ztest,nrow=1)
   
ntest <- dim(Ztest)[1]

## ESTIMATE THE DIRECTION BETA
########################################
BETA <- rep(0,length=r)
for (j in 1:ntrain)
{
   Wk <- WKernel[,j]
   
   fit <- wirrls(Y=Ytrain,Z=Ztrain,Lambda=Lambda,NbrIterMax=NbIterMax,WKernel=diag(c(Wk)))
   rm(Wk)
   
   #Check WIRRLS convergence
   if (fit$Cvg == 0)
   {
      return(list(Convergence=fit$Cvg))
   }
   
   cte <- fit$Coefficients[1,]
   gamma <- fit$Coefficients[-1,]
   
   #constante modification
   cte <- cte+Ztrain[j,-1]%*%gamma
   expcte <- exp(cte)
    
   gamma <- gamma*expcte/(1+expcte)-(gamma*expcte)*expcte/(1+expcte)^2
   BETA <- BETA + gamma/ntrain
}

rm(gamma)
## ESTIMATE hB
##########################################################
newsXtrain <- Ztrain[,-1]%*%BETA
if ((is.numeric(hB)==FALSE) || length(hB)!=1)
{
   res <-hplugin(newXtrain=newsXtrain,Ytrain=Ytrain,tau=0.1,intgsize=100,NbrIterMax=NbIterMax)
   hB <- res$h
}

## CLASSIFICATION STEP
##########################################################
newsXtest <- Ztest[,-1]%*%BETA 
cvg <- 1

Xikdes <- newsXtrain-c(newsXtest)
Wk <- exp((-Xikdes**2)/(2*hB^2))
Xikdes <-cbind(rep(1,length=ntrain),Xikdes)
    
#new fit over the FDR space
fit <- wirrls(Ytrain,Xikdes,Lambda=0,NbrIterMax=NbIterMax,WKernel=diag(c(Wk)))
cvg <- as.numeric((cvg==1) && (fit$Cvg==1))
ajEst = fit$Coefficients[1]


## PREDICTION STEP
##########################################################
piest <- 1/(1+exp(-ajEst))
hatY <- as.numeric(piest> 0.5)

## CONCLUDE
##############
return(list(Convergence=cvg,hatY=hatY))

}
