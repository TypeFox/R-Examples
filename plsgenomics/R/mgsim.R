### mgsim.R  (2006-01)
###
###                 GSIM for categorical data
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

mgsim <- function (Ytrain,Xtrain,Lambda,h,Xtest=NULL,NbIterMax=50)
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
##  h : real
##      indicates h parameter value
##  NbIterMax : positive integer
##      maximal number of iteration in the WIRRLS part


##    OUTPUT VARIABLES
##########################
## Ytest : vector ntest 
##      predicted label for ncomp
## beta : vector of length p*c (c=number of class -1)
##      estimated directions for the projection
## Coefficients : Coefficients corresponding to the step B of GSIM
## DeletedCol : vector
## if some covariables have nul variance, DeletedCol gives the 
## corresponding column number. Otherwise DeletedCol = NULL
## Cvg : 1 if convergence in MWIRRLS 0 otherwise

##  TEST ON INPUT VARIABLES
##############################
#On Xtrain
if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE)) {
 stop("Message from mgsim.R: Xtrain is not of valid type")}

if (dim(Xtrain)[2]==1) {
 stop("Message from mgsim.R: p=1 is not valid")}

ntrain <- dim(Xtrain)[1]
p <- dim(Xtrain)[2]

#On Xtest
if (is.null(Xtest)==FALSE) {

 if (is.vector(Xtest)==TRUE)
  {Xtest <- matrix(Xtest,nrow=1)}

 if ((is.matrix(Xtest)==FALSE)||(is.numeric(Xtest)==FALSE)) {
  stop("Message from mgsim.R: Xtest is not of valid type")}

 if (dim(Xtrain)[2]!=dim(Xtest)[2]) {
  stop("Message from mgsim.R: columns of Xtest and columns of Xtrain must be equal")}
  
 ntest <- dim(Xtest)[1] 
}

#On Ytrain
if ((is.vector(Ytrain)==FALSE)||(is.numeric(Ytrain)==FALSE)) {
 stop("Message from mgsim.R: Ytrain is not of valid type")}

if (length(Ytrain)!=ntrain) {
 stop("Message from mgsim.R: the length of Ytrain is not equal to the Xtrain row number")}

Ytrain <- Ytrain-1

if ((sum(floor(Ytrain)-Ytrain)!=0)||(sum(Ytrain<0)>0)){
 stop("Message from mgsim.R: Ytrain is not of valid type")}
 
c <- max(Ytrain)
eff<-rep(0,(c+1))
for (i in 0:c) {
    eff[(i+1)]<-sum(Ytrain==i)}
if (sum(eff==0)>0) {
 stop("Message from mgsim.R: there are empty classes")}

if (c==1) {
 stop("Message from mgsim.R: Ytrain is a binary vector, use gsim.R")}
 
#On hyper parameters 

if (is.vector(Lambda)==FALSE)
   stop("Message from mgsim.R: Lambda is not of valid type")
if (length(Lambda)!=1)
   stop("Message from mgsim.R: only one value can be specified for Lambda")
if ((is.numeric(Lambda)==FALSE)||(Lambda<0)){
 stop("Message from mgsim.R: Lambda is not of valid type")}

if (is.vector(h)==FALSE)
   stop("Message from mgsim.R: h is not of valid type")
if (length(h)!=1)
   stop("Message from mgsim.R: only one value can be specified for h")
if ((is.numeric(h)==FALSE)||(h<=0)){
 stop("Message from mgsim.R: h is not of valid type")}

if (is.vector(NbIterMax)==FALSE)
   stop("Message from mgsim.R: NbIterMax is not of valid type")
if (length(NbIterMax)!=1)
   stop("Message from mgsim.R: only one value can be specified for NbIterMax")
if ((is.numeric(NbIterMax)==FALSE)||(round(NbIterMax)-NbIterMax!=0)||(NbIterMax<1)){
 stop("Message from mgsim.R: NbIterMax is not of valid type")}

#Some initializations
r <- min(p,ntrain)
DeletedCol <- NULL
Cvg <- 1

##  MOVE IN THE REDUCED SPACE
################################
#  Standardize the Xtrain matrix
Sigma2train <- apply(Xtrain,2,var)*(ntrain-1)/ntrain
if (sum(Sigma2train==0)!=0){
    if (sum(Sigma2train==0)>(p-2)){
        stop("Message from mgsim.R: the procedure stops because number of predictor variables with no null variance is less than 1.")}
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


## ESTIMATE THE DIRECTIONS BETA
########################################
BETA <- matrix(0,nrow=r,ncol=c)
for (j in 1:ntrain) {
# Compute the Kernel matrix
 WK <- apply(exp(-(sXtrain-rep(1,ntrain)%*%t(sXtrain[j,]))^2/(2*h^2)),1,prod)
 WK <- WK/sum(WK)
 WKernel <- rep(0,ntrain*c)
 for (k in 1:ntrain) {
  WKernel[(1:c)+(k-1)*c] <- rep(WK[k],c)  
 }
 fit <- mwirrls(Y=Ytrain,Z=Zbloc,Lambda=Lambda,NbrIterMax=NbIterMax,WKernel=diag(c(WKernel)))
 rm(WKernel)
 rm(WK)
 #Check WIRRLS convergence
 if (fit$Cvg==0){
  warning("Message from mgsim : the algorithm did not converge in step A")
  Cvg <- 0}
 GAMMA <- matrix(fit$Coefficients,nrow=(r+1))
 cte <- GAMMA[1,]
 GAMMA <- GAMMA[-1,]
 cte=cte +sXtrain[j,]%*%GAMMA   
 cte <- exp(cte)
 GAMMA <- (GAMMA%*%diag(c(cte))-(GAMMA%*%t(cte))%*%cte/(1+sum(cte)))/(1+sum(cte))
 BETA <- BETA + GAMMA/ntrain
}
BETA <- sweep(BETA,2, sqrt(apply(BETA^2,2,sum)),FUN="/")
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
#Check WIRRLS convergence
if (fit$Cvg==0){
  warning("Message from mgsim : the algorithm did not converge in the second step (i.e. after projection)")
  Cvg <- 0}
  
## PREDICTION STEP
####################
hatY <- NULL
if (is.null(Xtest)==FALSE) {
 Zbloc <- Ztestbloc%*%B
 hatY <- apply(cbind(rep(0,ntest),matrix(Zbloc%*%fit$Coefficients,nrow=ntest,byrow=TRUE)),1,which.max)-1}



## CONCLUDE
##############

##Compute the estimated directions in the initial space

if (p>ntrain)
{beta <- diag(c(1/sqrt(Sigma2train)))%*%V%*%BETA}
if (p<=ntrain)
{beta <- diag(c(1/sqrt(Sigma2train)))%*%BETA}
List <- list(Ytest=(hatY+1),beta=beta,Coefficients=matrix(fit$Coefficients,nrow=2),DeletedCol=DeletedCol,Cvg=Cvg)
return(List)

}
