### mgsimaux.R  (2006-01)
###
###                 GSIM for categorical data
###              (procedure used in CV function)
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

mgsimaux <- function (Ytrain,Zbloc,sXtrain,Lambda,WKernel,Ztestbloc,NbIterMax=15) 
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
##  WKernel : matrix ntrain x ntrain
##      bandwith matrix
##  NbIterMax : positive integer
##      maximal number of iteration in the WIRRLS part


##    OUTPUT VARIABLES
##########################
## hatY : vector ntest 
##      predicted label for ncomp
## Cvg : 1 if convergence in WIRRLS 0 otherwise

c <- max(Ytrain)
r <- dim(Zbloc)[2]/c-1
ntrain <- dim(Zbloc)[1]/c
if (is.vector(Ztestbloc)==TRUE)
{Ztestbloc <- matrix(Ztestbloc,nrow=1)}
ntest <- dim(Ztestbloc)[1]/c


## ESTIMATE THE DIRECTIONS BETA
########################################
BETA <- matrix(0,nrow=r,ncol=c)
for (j in 1:ntrain) {
# Compute the Kernel matrix
 WKer <- rep(0,ntrain*c)
 for (k in 1:ntrain) {
  WKer[(1:c)+(k-1)*c] <- rep(WKernel[k,j],c)  
 }
 fit <- mwirrls(Y=Ytrain,Z=Zbloc,Lambda=Lambda,NbrIterMax=NbIterMax,WKernel=diag(c(WKer)))
 rm(WKer)
 #Check WIRRLS convergence
 if (fit$Cvg==0){
  return(list(Convergence=fit$Cvg))
  }
  
 GAMMA <- matrix(fit$Coefficients,nrow=(r+1))
 cte <- GAMMA[1,]
 GAMMA <- GAMMA[-1,]
 cte=cte +sXtrain[j,]%*%GAMMA   
 cte <- exp(cte)
 GAMMA <- (GAMMA%*%diag(c(cte))-(GAMMA%*%t(cte))%*%cte/(1+sum(cte)))/(1+sum(cte))
 BETA <- BETA + GAMMA/ntrain
}

rm(GAMMA)

## ESTIMATE COEFFICIENTS FOR THE DESIGN AFTER PROJECTING 
##########################################################
B <- matrix(0,nrow=(r+1)*c,ncol=2*c)
for (k in 1:c) {
 B[((k-1)*(r+1)+1),(2*(k-1)+1)] <- 1
 B[((k-1)*(r+1)+2):(k*(r+1)),2*k] <- BETA[,k]
}
Zbloc <- Zbloc%*%B
fit <- mwirrls(Y=Ytrain,Z=Zbloc,Lambda=0,NbrIterMax=NbIterMax,WKernel=diag(rep(1,ntrain*c)))
#Check WIRRLS convergence ? faire un message uniquement pour un retour 11111???  
 if (fit$Cvg==0){
  return(list(Convergence=fit$Cvg))
  }

## PREDICTION STEP
####################
Zbloc <- Ztestbloc%*%B
hatY <- apply(cbind(rep(0,ntest),matrix(Zbloc%*%fit$Coefficients,nrow=ntest,byrow=TRUE)),1,which.max)-1



## CONCLUDE
##############
return(list(Convergence=fit$Cvg,hatY=hatY))

}
