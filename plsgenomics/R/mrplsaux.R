### mrplsaux.R  (2006-01)
###
###    Ridge Partial Least square for categorical data
###          (procedure used in CV function)
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

mrplsaux <- function (Ytrain,Zbloc,Lambda,ncomp,Ztestbloc,NbIterMax=15)
{

##    INPUT VARIABLES
#########################
##  Zbloc   : matrix ntrain*c x (r+1)*c
##      train design matrix
##  Ytrain   : vector ntrain
##      response variable {0,1,...,c}-valued vector
##  Ztestbloc   : matrix ntest*c x (r+1)*c
##      test design matrix
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
## for ncomp=k for ncomp !=0,1
## Cvg : 1 if convergence in WIRRLS 0 otherwise

c <- max(Ytrain)
r <- dim(Zbloc)[2]/c-1
ntrain <- dim(Zbloc)[1]/c
if (is.vector(Ztestbloc)==TRUE)
{Ztestbloc <- matrix(Ztestbloc,nrow=1)}
ntest <- dim(Ztestbloc)[1]/c


## RUN RPLS IN THE REDUCED SPACE
########################################

fit <- mwirrls(Y=Ytrain,Z=Zbloc,Lambda=Lambda,NbrIterMax=NbIterMax,WKernel=diag(rep(1,ntrain*c)))
#Check WIRRLS convergence
if (fit$Cvg==0)
stop("Message from rpls : WIRRLS did not converge; try another Lambda value")

if (ncomp==0) #Ridge procedure
{GAMMA <- fit$Coefficients
 hatY <- apply(cbind(rep(0,ntest),matrix(Ztestbloc%*%GAMMA,nrow=ntest,byrow=TRUE)),1,which.max)-1}

if (ncomp!=0) { 
#Compute Weight and pseudo variable
#Pseudovar = Eta + W^-1 Psi
Eta <- Zbloc%*%fit$Coefficients
mu <- rep(0,length(Eta))
W <- matrix(0,length(mu),length(mu))
for (kk in 1:ntrain) {
 mu[c*(kk-1)+(1:c)] <- exp(Eta[c*(kk-1)+(1:c)])/(1+sum(exp(Eta[c*(kk-1)+(1:c)])))
 Blocmu <- mu[c*(kk-1)+(1:c)]
 BlocW <- -Blocmu%*%t(Blocmu)
 BlocW <- BlocW+diag(Blocmu)
 W[c*(kk-1)+(1:c),c*(kk-1)+(1:c)] <- BlocW
 }
Psi <- fit$Ybloc-mu


## Run PLS

# W-Center the sXtrain and pseudo variable
col <- seq(from=1, to=c*(r+1), by= (r+1))
Xbloc <- Zbloc[,-col]
Cte <- Zbloc[,col]  
# Weighted centering of Pseudo variable
H <- t(Cte)%*%W%*%Cte
WMeanPseudoVar <- solve(H,t(Cte)%*%(W%*%Eta+Psi))
WCtrPsi <- Psi
WCtrEta <- Eta-Cte%*%WMeanPseudoVar
# Weighted centering of sXtrain
WMeansXtrain <- solve(H,t(Cte)%*%W%*%Xbloc)
WCtrsXtrain <- Xbloc-Cte%*%WMeansXtrain
rm(H)

#Initialize some variables
PsiAux <- diag(c(rep(1,r*c)))
E <- WCtrsXtrain
f1 <- WCtrEta
f2 <- WCtrPsi
Omega <- matrix(0,nrow=r*c,ncol=ncomp)
Scores <- matrix(0,nrow=ntrain*c,ncol=ncomp)
TildePsi <- matrix(0,nrow=r*c,ncol=ncomp)
Loadings <- matrix(0,nrow=r*c,ncol=ncomp)
qcoeff <- vector(ncomp,mode="numeric")
GAMMA <- matrix(0,nrow=c*(r+1),ncol=ncomp)
hatY <- matrix(0,nrow=ntest,ncol=ncomp)

#WPLS loop
for (count in 1:ncomp)
{Omega[,count]<-t(E)%*%(W%*%f1+f2)
 #Score vector
 t<-E%*%Omega[,count]
 c1<-t(Omega[,count])%*%t(E)%*%W%*%E%*%Omega[,count]
 Scores[,count]<-t
 TildePsi[,count] <- PsiAux%*%Omega[,count]
 #Deflation of X
 Loadings[,count]<-t(t(t)%*%W%*%E)/c1[1,1]
 E<-E-t%*%t(Loadings[,count])
 #Deflation of f1
 qcoeff[count]<-t(W%*%f1+f2)%*%t/c1[1,1]
 f1 <- f1-qcoeff[count]*t
 #Recursive definition of RMatrix
 PsiAux<-PsiAux%*%(diag(c(rep(1,c*r)))-Omega[,count]%*%t(Loadings[,count]))
 #Express regression coefficients w.r.t. the columns of [1 sX] for ncomp=count
 if (count==1)
  {GAMMA[-col,count]<-TildePsi[,1:count]%*%t(c(qcoeff[1:count]))}
 if (count!=1)
  {GAMMA[-col,count]<-TildePsi[,1:count]%*%qcoeff[1:count]}
 GAMMA[col,count]=WMeanPseudoVar-WMeansXtrain%*%GAMMA[-col,count]
 #classification step
 ETA <- cbind(rep(0,ntest),matrix(Ztestbloc%*%GAMMA[,count],nrow=ntest,byrow=TRUE))
 hatY[,count] <- apply(ETA,1,which.max)-1}

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
