### mrpls.R  (2006-01)
###
###    Ridge Partial Least square for categorical data
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

mrpls <- function (Ytrain,Xtrain,Lambda,ncomp,Xtest=NULL,NbIterMax=50)
{

##    INPUT VARIABLES
#########################
##  Xtrain   : matrix ntrain x p
##      train data matrix
##  Ytrain   : vector ntrain
##      response variable {1,...,c+1}-valued vector
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
## Ytest : vector ntest 
##      predicted label for ncomp
## Coefficients : matrix p+1 x c 
##                regression coefficients w.r.t. the columns of [1 Xtest]
## DeletedCol : vector
## if some covariables have nul variance, DeletedCol gives the 
## corresponding column number. Otherwise DeletedCol = NULL

##  TEST ON INPUT VARIABLES
##############################
#On Xtrain
if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE)) {
 stop("Message from mrpls.R: Xtrain is not of valid type")}

if (dim(Xtrain)[2]==1) {
 stop("Message from mrpls.R: p=1 is not valid")}

ntrain <- dim(Xtrain)[1]
p <- dim(Xtrain)[2]

#On Xtest
if (is.null(Xtest)==FALSE) {

 if (is.vector(Xtest)==TRUE)
  {Xtest <- matrix(Xtest,nrow=1)}

 if ((is.matrix(Xtest)==FALSE)||(is.numeric(Xtest)==FALSE)) {
  stop("Message from mrpls.R: Xtest is not of valid type")}

 if (dim(Xtrain)[2]!=dim(Xtest)[2]) {
  stop("Message from mrpls.R: columns of Xtest and columns of Xtrain must be equal")}
  
 ntest <- dim(Xtest)[1] 
}

#On Ytrain
if ((is.vector(Ytrain)==FALSE)||(is.numeric(Ytrain)==FALSE)) {
 stop("Message from mrpls.R: Ytrain is not of valid type")}

if (length(Ytrain)!=ntrain) {
 stop("Message from mrpls.R: the length of Ytrain is not equal to the Xtrain row number")}

Ytrain <- Ytrain-1

if ((sum(floor(Ytrain)-Ytrain)!=0)||(sum(Ytrain<0)>0)){
 stop("Message from mrpls.R: Ytrain is not of valid type")}
 
c <- max(Ytrain)
eff<-rep(0,(c+1))
for (i in 0:c) {
    eff[(i+1)]<-sum(Ytrain==i)}
if (sum(eff==0)>0) {
 stop("Message from mrpls.R: there are empty classes")}

if (c==1) {
 stop("Message from mrpls.R: Ytrain is a binary vector, use rpls.R")}
 
#On hyper parameters 

if ((is.numeric(Lambda)==FALSE)||(Lambda<0)){
 stop("Message from mrpls.R: Lambda is not of valid type")}

if ((is.numeric(ncomp)==FALSE)||(round(ncomp)-ncomp!=0)||(ncomp<0)){
 stop("Message from mrpls.R: ncomp is not of valid type")}

if ((is.numeric(NbIterMax)==FALSE)||(round(NbIterMax)-NbIterMax!=0)||(NbIterMax<1)){
 stop("Message from mrpls.R: NbIterMax is not of valid type")}

 
#Some initializations
r <- min(p,ntrain)
c <- max(Ytrain)
DeletedCol <- NULL


##  MOVE IN THE REDUCED SPACE
################################
#  Standardize the Xtrain matrix
Sigma2train <- apply(Xtrain,2,var)*(ntrain-1)/ntrain
if (sum(Sigma2train==0)!=0){
    if (sum(Sigma2train==0)>(p-2)){
        stop("Message from mrpls.R: the procedure stops because number of predictor variables with no null variance is less than 1.")}
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

rm(Xtrain)
#Compute Zblock 
Z <- cbind(rep(1,ntrain),sXtrain)
Zbloc <- matrix(0,nrow=ntrain*c,ncol=c*(r+1))
 
if (is.null(Xtest)==FALSE) {
 sXtest <- sweep(Xtest,2,MeanXtrain,FUN="-")
 sXtest <- sweep(sXtest,2,sqrt(Sigma2train),FUN="/")
 if (p>ntrain)
  {sXtest <- sXtest%*%V}
 Xtest  <- 0
 Zt <- cbind(rep(1,ntest),sXtest)
 Ztestbloc <- matrix(0,nrow=ntest*c,ncol=c*(r+1))}
 
for (cc in 1:c) {
  row <- (0:(ntrain-1))*c+cc
  col <- (r+1)*(cc-1)+1:(r+1)
  Zbloc[row,col] <- Z
  if (is.null(Xtest)==FALSE) {
    row <- (0:(ntest-1))*c+cc
    Ztestbloc[row,col] <- Zt}
}       
rm(Z)
Zt <- NULL

## RUN RPLS IN THE REDUCED SPACE
########################################

fit <- mwirrls(Y=Ytrain,Z=Zbloc,Lambda=Lambda,NbrIterMax=NbIterMax,WKernel=diag(rep(1,ntrain*c)))
#Check WIRRLS convergence
if (fit$Cvg==0)
stop("Message from rpls : WIRRLS did not converge; try another Lambda value")

if (ncomp==0) #Ridge procedure
{GAMMA <- fit$Coefficients
 if (is.null(Xtest)==FALSE) {
 hatY <- apply(cbind(rep(0,ntest),matrix(Ztestbloc%*%GAMMA,nrow=ntest,byrow=TRUE)),1,which.max)-1}
}

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
index <- 1:(c*(r+1))
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
Omega <- matrix(0,r*c,ncomp)
Scores <- matrix(0,ntrain*c,ncomp)
TildePsi <- matrix(0,r*c,ncomp)
Loadings <- matrix(0,r*c,ncomp)
qcoeff <- vector(ncomp,mode="numeric")
GAMMA <- matrix(0,nrow=c*(r+1),ncol=ncomp)
if (is.null(Xtest)==FALSE) {
 hatY <- matrix(0,nrow=ntest,ncol=ncomp)}

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
 if (is.null(Xtest)==FALSE) {
   hatY[,count] <- apply(cbind(rep(0,ntest),matrix(Ztestbloc%*%GAMMA[,count],nrow=ntest,byrow=TRUE)),1,which.max)-1}
 }  

}


## CONCLUDE
##############

##Compute the coefficients w.r.t. [1 X]
if (ncomp!=0)
{GAMMA <- GAMMA[,ncomp]}
GAMMA <- t(matrix(GAMMA,nrow=c,byrow=TRUE))
Coefficients <- t(matrix(0,nrow=c,ncol=(p+1)))
if (p>ntrain)
{Coefficients[-1,] <- diag(c(1/sqrt(Sigma2train)))%*%V%*%GAMMA[-1,]}
if (p<=ntrain)
{Coefficients[-1,] <- diag(c(1/sqrt(Sigma2train)))%*%GAMMA[-1,]}
Coefficients[1,] <- GAMMA[1,]-MeanXtrain%*%Coefficients[-1,]
List <- list(Coefficients=Coefficients,Ytest=NULL,DeletedCol=DeletedCol)
if (is.null(Xtest)==FALSE) {
    if ((ncomp==0)|(ncomp==1))
        {List <- list(Coefficients=Coefficients,Ytest=(hatY+1),DeletedCol=DeletedCol)}
    if ((ncomp!=0)&(ncomp!=1))
        {colnames(hatY)=1:ncomp
         rownames(hatY)=1:ntest
         List <- list(Coefficients=Coefficients,hatY=(hatY+1),Ytest=(hatY[,ncomp]+1),DeletedCol=DeletedCol)}
}
return(List)

}

