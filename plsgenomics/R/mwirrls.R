### mwirrls.R  (2006-01)
###
###    Local Weighted Iteratively Reweighted Ridge Least Squares (IRRLS)
###                      for categorical data
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

mwirrls <- function(Y,Z,Lambda=0,NbrIterMax=15,Threshold=10^(-12),WKernel) 
{
##    IN     
##########
##      c=NbrClass-1
##  Z : Design matrix (cn) x c(p+1)
##          block data matrix
##  Y : Response Vector
##          vector ntrain 
##  Lambda :  coefficient of the ridge penalty
##          real
##  NbrIterMax : Maximal number of iterations
##          positive integer
##  Threshold : Used for the stopping rule.
##          real
##  WKernel : Kernel Weigth
##          matrix cn x cn 

##    OUT     
############
##  out :  structure that contains the fields
##      Gamma : vector of the regression coefficients w.r.t the design
##      matrix.
##      Cvg : Cvg=1 if the algorithm has converged otherwise 0.
##      Ybloc : block Response Vector ntrain*c

c <- max(Y)
dZ2 <- dim(Z)[2]
n <- length(Y)
p <- dZ2/c-1

Ybloc <- rep(0,(n*c));  
for (cc in 1:c)
{ff <- which(Y==cc)  
 Ybloc[(ff-1)*c+cc]=rep(1,length(ff))}

R <- matrix(0,dZ2,dZ2) 
diag(R) <- rep(1,dZ2)
grid <- seq(from=1, to=dZ2, by= (p+1))
R[grid,grid]<-0 



##  THE NR ALGORITHM
##################################
##  1. Initialize the parameter

mu <- (1+2*Ybloc)/(c+3)
Eta <- rep(0,length(mu))
WNR <- matrix(0,length(mu),length(mu))

for (kk in 1:n) {
 ss <- 1-sum(mu[c*(kk-1)+(1:c)])
 Blocmu <- mu[c*(kk-1)+(1:c)]
 BlocW <- -Blocmu%*%t(Blocmu)
 BlocW <- BlocW+diag(Blocmu)
 WNR[c*(kk-1)+(1:c),c*(kk-1)+(1:c)] <- BlocW
 Eta[c*(kk-1)+(1:c)] <- log(Blocmu)-log(ss)
}
             
WT <- WKernel%*%WNR           
H <- t(Z)%*%WT%*%Z+Lambda*R    
trysolve<-try(solve(H,t(Z)%*%WKernel%*%(WNR%*%Eta+(Ybloc-mu))), silent = TRUE)
if (sum(is.nan(trysolve))>0) {
 trysolve <- NULL}

if (is.matrix(trysolve)==FALSE) {
  Gamma <- rep(1,dZ2) 
  StopRule <- 1   
  NbrIter <- NbrIterMax+1
  Illcond <- 1
  Separation <- 0
  Cvg <- 0}
if (is.matrix(trysolve)==TRUE)  {
  StopRule <- 0 
  NbrIter <- 0
  Illcond <- 0
  Separation <- 0}              
            
##  2. Newton-Raphson loop
 
while (StopRule==0)
{#Increment the iterations number 
 NbrIter <- NbrIter + 1
 #Udapte Gamma
 if (NbrIter==1)
  {Gamma <- trysolve}
 if (NbrIter>1)
  {Gamma <- Gamma+trysolve}           
 #Udapte Eta            
 Eta <- Z%*%Gamma  
 Eta[Eta>700] <- rep(700,sum(Eta>700)) 
 #Udapte mu and WNR
 for (kk in 1:n) {
 mu[c*(kk-1)+(1:c)] <- exp(Eta[c*(kk-1)+(1:c)])/(1+sum(exp(Eta[c*(kk-1)+(1:c)])))
 Blocmu <- mu[c*(kk-1)+(1:c)]
 BlocW <- -Blocmu%*%t(Blocmu)
 BlocW <- BlocW+diag(Blocmu)
 WNR[c*(kk-1)+(1:c),c*(kk-1)+(1:c)] <- BlocW
 }
 #Udapte total Weight   
 WT <- WNR%*%WKernel
 #Udapte Gradient
 Gradient <- t(Z)%*%WKernel%*%(Ybloc-mu)-Lambda*R%*%Gamma 
 #Udapte H           
 H <- t(Z)%*%WT%*%Z+Lambda*R   
 trysolve <- try(solve(H,Gradient), silent = TRUE)      
 if (sum(is.nan(trysolve))>0) {
 trysolve <- NULL}
        
 #Compute the StopRule
 #on the (quasi)-separation detection (for Lambda=0)
 if (Lambda==0)
  {Separation <- as.numeric(sum((apply(cbind(rep(0,n),matrix(Eta,nrow=n,byrow=TRUE)),1,which.max)-1)==Y)==n)}
 #on the conditioning of matrix
 Illcond <- as.numeric(is.matrix(trysolve)==FALSE)
 #on the convergence
 Cvg <- as.numeric(sqrt(sum(abs(Gradient)^2))<=Threshold)
 StopRule <- max(Cvg,as.numeric(NbrIter>=NbrIterMax),Separation,Illcond) 
}

##  CONCLUSION
###############        

list <- list(Coefficients=Gamma,Cvg=max(Cvg,Separation),Ybloc=Ybloc)
return(list)
}
