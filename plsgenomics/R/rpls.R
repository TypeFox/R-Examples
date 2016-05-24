### rpls.R  (2006-01)
###
###    Ridge Partial Least square for binary data
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

rpls <- function (Ytrain,Xtrain,Lambda,ncomp,Xtest=NULL,NbIterMax=50)
{

##    INPUT VARIABLES
#########################
##  Xtrain   : matrix ntrain x p
##      train data matrix
##  Ytrain   : vector ntrain
##      response variable {1,2}-valued vector
##  Xtest   : NULL or matrix ntest x p 
##      if no NULL Xtest is the test data matrix
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
## Ytest : matrix ntest x 1
##      predicted label for ncomp
## Coefficients : vector p+1 x 1
##                regression coefficients w.r.t. the columns of [1 Xtest]
## DeletedCol : vector
## if some covariables have nul variance, DeletedCol gives the 
## corresponding column number. Otherwise DeletedCol = NULL

##  TEST ON INPUT VARIABLES
##############################
#On Xtrain
if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE)) {
 stop("Message from rpls.R: Xtrain is not of valid type")}

if (dim(Xtrain)[2]==1) {
 stop("Message from rpls.R: p=1 is not valid")}

ntrain <- dim(Xtrain)[1]
p <- dim(Xtrain)[2]

#On Xtest
if (is.null(Xtest)==FALSE) {

 if (is.vector(Xtest)==TRUE)
  {Xtest <- matrix(Xtest,nrow=1)}

 if ((is.matrix(Xtest)==FALSE)||(is.numeric(Xtest)==FALSE)) {
  stop("Message from rpls.R: Xtest is not of valid type")}

 if (dim(Xtrain)[2]!=dim(Xtest)[2]) {
  stop("Message from rpls.R: columns of Xtest and columns of Xtrain must be equal")}
  
 ntest <- dim(Xtest)[1] 
}

#On Ytrain
if ((is.vector(Ytrain)==FALSE)||(is.numeric(Ytrain)==FALSE)) {
 stop("Message from rpls.R: Ytrain is not of valid type")}

if (length(Ytrain)!=ntrain) {
 stop("Message from rpls.R: the length of Ytrain is not equal to the Xtrain row number")}

Ytrain <- Ytrain-1

if ((sum(floor(Ytrain)-Ytrain)!=0)||(sum(Ytrain<0)>0)){
 stop("Message from rpls.R: Ytrain is not of valid type")}
 
c <- max(Ytrain)
if (c!=1) { 
 stop("Message from rpls.R: Ytrain is not of valid type")}

eff<-rep(0,2)
for (i in 0:1) {
    eff[i+1]<-sum(Ytrain==i)}
if (sum(eff==0)>0) {
 stop("Message from rpls.R: there are empty classes")}

 
#On hyper parameters 

if ((is.numeric(Lambda)==FALSE)||(Lambda<0)){
 stop("Message from rpls.R: Lambda is not of valid type")}

if ((is.numeric(ncomp)==FALSE)||(round(ncomp)-ncomp!=0)||(ncomp<0)){
 stop("Message from rpls.R: ncomp is not of valid type")}

if ((is.numeric(NbIterMax)==FALSE)||(round(NbIterMax)-NbIterMax!=0)||(NbIterMax<1)){
 stop("Message from rpls.R: NbIterMax is not of valid type")}

 
#Some initializations
r <- min(p,ntrain)
DeletedCol <- NULL

##  MOVE IN THE REDUCED SPACE
################################
#  Standardize the Xtrain matrix

Sigma2train <- apply(Xtrain,2,var)*(ntrain-1)/ntrain
if (sum(Sigma2train==0)!=0){
    if (sum(Sigma2train==0)>(p-2)){
        stop("Message from rpls.R: the procedure stops because number of predictor variables with no null variance is less than 1.")}
    warning("There are covariables with nul variance")
    Xtrain <- Xtrain[,which(Sigma2train!=0)]
    Xtest <- Xtest[,which(Sigma2train!=0)]
    if (is.vector(Xtest)==TRUE)
     {Xtest <- matrix(Xtest,nrow=1)}
    index <- 1:p
    DeletedCol <- index[which(Sigma2train==0)]
    Sigma2train <-Sigma2train[which(Sigma2train!=0)]
    p <- dim(Xtrain)[2]
    r <- min(p,ntrain)}
MeanXtrain <- apply(Xtrain,2,mean)
sXtrain <- sweep(Xtrain,2,MeanXtrain,FUN="-")
sXtrain <- sweep(sXtrain,2,sqrt(Sigma2train),FUN="/")

#Compute the svd when necessary
if (p>ntrain)
{svd.sXtrain <- svd(t(sXtrain))
 r<-length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
 V <- svd.sXtrain$u[,1:r]
 D <- diag(c(svd.sXtrain$d[1:r]))
 U <- svd.sXtrain$v[,1:r]
 sXtrain <- U%*%D
 rm(D)
 rm(U)
 rm(svd.sXtrain)}
 
if (is.null(Xtest)==FALSE) {
 sXtest <- sweep(Xtest,2,MeanXtrain,FUN="-")
 sXtest <- sweep(sXtest,2,sqrt(Sigma2train),FUN="/")
 if (p>ntrain)
  {sXtest <- sXtest%*%V}
 Xtest <- 0}

rm(Xtrain)


## RUN RPLS IN THE REDUCED SPACE
########################################

fit <- wirrls(Y=Ytrain,Z=cbind(rep(1,ntrain),sXtrain),Lambda=Lambda,NbrIterMax=NbIterMax,WKernel=diag(rep(1,ntrain)))
#Check WIRRLS convergence
if (fit$Cvg==0)
stop("Message from rpls : WIRRLS did not converge; try another Lambda value")

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
  {GAMMA[-1,count]<-TildePsi[,1:count]%*%t(c(qcoeff[1:count]))}
 if (count!=1)
  {GAMMA[-1,count]<-TildePsi[,1:count]%*%qcoeff[1:count]}
 GAMMA[1,count]=WMeanPseudoVar-WMeansXtrain%*%GAMMA[-1,count]}}

## CLASSIFICATION STEP
#######################
if (is.null(Xtest)==FALSE) {
 hatY <- cbind(rep(1,ntest),sXtest)%*%GAMMA
 hatY <- (hatY>0)+0}


## CONCLUDE
##############


##Compute the coefficients w.r.t. [1 X]
if (ncomp!=0)
{GAMMA <- GAMMA[,ncomp]}
Coefficients <- rep(0,p+1)
if (p>ntrain)
{Coefficients[-1] <- diag(c(1/sqrt(Sigma2train)))%*%V%*%GAMMA[-1]}
if (p<=ntrain)
{Coefficients[-1] <- diag(c(1/sqrt(Sigma2train)))%*%GAMMA[-1]}
Coefficients[1] <- GAMMA[1]-MeanXtrain%*%Coefficients[-1]
List <- list(Coefficients=Coefficients,Ytest=NULL,DeletedCol=DeletedCol)
if (is.null(Xtest)==FALSE) {
    if ((ncomp==0)|(ncomp==1))
        {List <- list(Coefficients=Coefficients,Ytest=(hatY[,1]+1),DeletedCol=DeletedCol)}
    if ((ncomp!=0)&(ncomp!=1))
        {colnames(hatY)=1:ncomp
         rownames(hatY)=1:ntest
         List <- list(Coefficients=Coefficients,hatY=(hatY+1),Ytest=(hatY[,ncomp]+1),DeletedCol=DeletedCol)}
}
return(List)

}
