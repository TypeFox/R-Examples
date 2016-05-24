### rplsaux.R  (2006-01)
###
###    Ridge Partial Least square for binary data
###       (procedure used in CV function)
###
### Copyright 2006-01 Sophie Lambert-Lacroix
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

rplsaux <- function (Ytrain,sXtrain,Lambda,ncomp,sXtest,NbIterMax=50)
{

##    INPUT VARIABLES
#########################
##  sXtrain   : matrix ntrain x r
##      train data matrix
##  Ytrain   : vector ntrain
##      response variable {0,1}-valued vector
##  sXtest   : matrix ntest x r
##      test data matrix
##  Lambda : real
##      value for the regularization parameter Lambda
##  NbIterMax : positive integer
##      maximal number of iteration in the WIRRLS part
##  ncomp : maximal number of PLS components
##          0 = Ridge


##    OUTPUT VARIABLES
##########################
## hatY : matrix of size ntest x ncomp in such a way
## that the kieme column corresponds to the result
## for ncomp=k
## Cvg : 1 if convergence in WIRRLS 0 otherwise

r <- dim(sXtrain)[2]
ntrain <- dim(sXtrain)[1]
if (is.vector(sXtest)==TRUE)
{sXtest <- matrix(sXtest,nrow=1)}
ntest <- dim(sXtest)[1]

## RUN RPLS IN THE REDUCED SPACE
########################################

fit <- wirrls(Y=Ytrain,Z=cbind(rep(1,ntrain),sXtrain),Lambda=Lambda,NbrIterMax=NbIterMax,WKernel=diag(rep(1,ntrain)))

if (fit$Cvg==1) {

if (ncomp==0) #Ridge procedure
{GAMMA <- fit$Coefficients}

if (ncomp!=0) {
#Compute Weight and pseudo variable
#Pseudovar = Eta + W^-1 Psi
Eta <- cbind(rep(1,ntrain),sXtrain)%*%fit$Coefficients
mu<-1/(1+exp(-Eta))
diagW <- mu*(1-mu)
W <- diag(c(diagW))
Psi <- Ytrain-mu

## Run PLS
# W-Center the sXtrain and pseudo variable
Sum=sum(diagW)
# Weighted centering of Pseudo variable
WMeanPseudoVar <- sum(W%*%Eta+Psi)/Sum
WCtrPsi <- Psi
WCtrEta <- Eta-c(WMeanPseudoVar)
# Weighted centering of sXtrain
WMeansXtrain <- t(diagW)%*%sXtrain/Sum
WCtrsXtrain <- sXtrain-rep(1,ntrain)%*%WMeansXtrain

#Initialize some variables
PsiAux <- diag(c(rep(1,r)))
E <- WCtrsXtrain
f1 <- WCtrEta
f2 <- WCtrPsi
Omega <- matrix(0,r,ncomp)
Scores <- matrix(0,ntrain,ncomp)
TildePsi <- matrix(0,r,ncomp)
Loadings <- matrix(0,r,ncomp)
qcoeff <- vector(ncomp,mode="numeric")
GAMMA <- matrix(0,nrow=(r+1),ncol=ncomp)

#WPLS loop
for (count in 1:ncomp)
{Omega[,count]<-t(E)%*%(W%*%f1+f2)
 #Score vector
 t<-E%*%Omega[,count]
 c<-t(Omega[,count])%*%t(E)%*%W%*%E%*%Omega[,count]
 Scores[,count]<-t
 TildePsi[,count] <- PsiAux%*%Omega[,count]
 #Deflation of X
 Loadings[,count]<-t(t(t)%*%W%*%E)/c[1,1]
 E<-E-t%*%t(Loadings[,count])
 #Deflation of f1
 qcoeff[count]<-t(W%*%f1+f2)%*%t/c[1,1]
 f1 <- f1-qcoeff[count]*t
 #Recursve definition of RMatrix
 PsiAux<-PsiAux%*%(diag(c(rep(1,r)))-Omega[,count]%*%t(Loadings[,count]))
 #Express regression coefficients w.r.t. the columns of [1 sX] for ncomp=count
 if (count==1)
  {GAMMA[-1,count] <- TildePsi[,1:count]%*%t(c(qcoeff[1:count]))}
 if (count!=1)
  {GAMMA[-1,count] <- TildePsi[,1:count]%*%qcoeff[1:count]}
 GAMMA[1,count] <- WMeanPseudoVar-WMeansXtrain%*%GAMMA[-1,count]}
 }

## CLASSIFICATION STEP
#######################

ETA <- cbind(rep(1,ntest),sXtest)%*%GAMMA
hatY <- (ETA>0)+0
}


## CONCLUDE
##############

Cvg <- fit$Cvg

if (Cvg==1)
{List <- list(Convergence=Cvg,hatY=hatY)}
if (Cvg==0)
{List <- list(Convergence=Cvg)}

return(List)

}
