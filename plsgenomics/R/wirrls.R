### wirrls.R  (2006-01)
###
###    Local Weighted Iteratively Reweighted Ridge Least Squares (IRRLS)
###                      for binary data
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

wirrls <- function(Y,Z,Lambda=0,NbrIterMax=15,Threshold=10^(-12),WKernel) 
{
##    IN     
##########
##  Z : Design matrix
##          matrix n x (p+1)
##  Y : vector n 
##      response variable {0,1}-valued vector
##  Lambda :  coefficient of the ridge penalty
##          real
##  NbrIterMax : Maximal number of iterations
##          integer
##  Threshold : Used for the stopping rule.
##          real
##  WKernel : Kernel Weigth
##          matrix n x n 

##    OUT     
############
##  out :  structure that contains the fields
##      Gamma : vector of the regression coefficients w.r.t the design
##      matrix Z=[1 X].
##      Cvg : Cvg=1 if the algorithm has converged otherwise 0.

dZ2 <- dim(Z)[2]
n <- dim(Z)[1]
p <- dZ2-1
R<-matrix(0,dZ2,dZ2) 
R[2:dZ2,2:dZ2]<-diag(1,nrow=p) 

##  THE NR ALGORITHM
##################################
##  1. Initialize the parameter

c <- log(3) 
Eta <- c*Y-c*(1-Y)                
mu <- 1/(1+exp(-Eta))             
DiagWNR <- mu*(1-mu)              
WT <- diag(c(DiagWNR))%*%WKernel          
H <- t(Z)%*%WT%*%Z+Lambda*R                
trysolve<-try(solve(H,t(Z)%*%WKernel%*%(diag(c(DiagWNR))%*%Eta+(Y-mu))), silent = TRUE)
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
 #Udapte mu
 mu<-1/(1+exp(-Eta))   
 #Udapte Weight   
 DiagWNR <- mu*(1-mu)  
 WT <- diag(c(DiagWNR))%*%WKernel
 #Udapte Gradient
 Gradient <- t(Z)%*%WKernel%*%(Y-mu)-Lambda*R%*%Gamma
 #Udapte H           
 H <- t(Z)%*%WT%*%Z+Lambda*R   
 trysolve <- try(solve(H,Gradient), silent = TRUE)  
 if (sum(is.nan(trysolve))>0) {
 trysolve <- NULL}       
          
 #Compute the StopRule
 #on the (quasi)-separation detection (for Lambda=0)
 if (Lambda==0)
  {Separation <- as.numeric(sum((Eta>0)+0==Y)==n)}
 #on the conditioning of matrix
 Illcond <- as.numeric(is.matrix(trysolve)==FALSE) 
 #on the convergence
 Cvg <- as.numeric(sqrt(sum(abs(Gradient)^2))<=Threshold)
 StopRule <- max(Cvg,as.numeric(NbrIter>=NbrIterMax),Separation,Illcond) 
}

##  CONCLUSION
###############        

list <- list(Coefficients=Gamma,Cvg=max(Cvg,Separation))
return(list)
}
