### mgsim.cv.R  (2006-01)
###
###    Determination by Cross-validation of MGSIM
###        hyper-parameters for categorical data
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

mgsim.cv <- function (Ytrain,Xtrain,LambdaRange,hRange,NbIterMax=50)
{

##    INPUT VARIABLES
#########################
##  Xtrain   : matrix ntrain x p
##      train data matrix
##  Ytrain   : vector ntrain
##      response variable {1,...c+1}-valued vector
##  LambdaRange : vector nLambda
##      possible values for the regularization parameter Lambda
##  NbIterMax : positive integer
##      maximal number of iteration in the MWIRRLS part
##  hRandge : vector nh
##      possible values for the bandwidth parameter


##    OUTPUT VARIABLES
##########################
## Lambda : optimal regularization parameter Lambda
## h : optimal bandwidth parameter


##  TEST ON INPUT VARIABLES
##############################
#On Xtrain
if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE)) {
 stop("Message from mgsim.cv.R: Xtrain is not of valid type")}

if (dim(Xtrain)[2]==1) {
 stop("Message from mgsim.cv.R: p=1 is not valid")}

ntrain <- dim(Xtrain)[1]

#On Ytrain
if ((is.vector(Ytrain)==FALSE)||(is.numeric(Ytrain)==FALSE)) {
 stop("Message from mgsim.cv.R: Ytrain is not of valid type")}

if (length(Ytrain)!=ntrain) {
 stop("Message from mgsim.cv.R: the length of Ytrain is not equal to the Xtrain row number")}

Ytrain <- Ytrain-1
if ((sum(floor(Ytrain)-Ytrain)!=0)||(sum(Ytrain<0)>0)){
 stop("Message from mgsim.cv.R: Ytrain is not of valid type")}

c <- max(Ytrain)
eff<-rep(0,(c+1))
for (i in 0:c) {
    eff[(i+1)]<-sum(Ytrain==i)}
if (sum(eff<=1)>0) {
 stop("Message from mgsim.cv.R: there are not enough samples for each class")}

if (c==1) {
 stop("Message from mgsim.cv.R: Ytrain is a binary vector, use gsim.cv.R")}

#On hyper parameters range

if ((is.numeric(LambdaRange)==FALSE)||(is.vector(LambdaRange)==FALSE)||(sum(LambdaRange<0)>0)){
 stop("Message from mgsim.cv.R: LambdaRange is not of valid type")}
 
if ((is.numeric(hRange)==FALSE)||(is.vector(hRange)==FALSE)||(sum(hRange<=0)>0)){
 stop("Message from mgsim.cv.R: hRange is not of valid type")}

if ((is.numeric(NbIterMax)==FALSE)||(round(NbIterMax)-NbIterMax!=0)||(NbIterMax<1)){
 stop("Message from mgsim.cv.R: NbIterMax is not of valid type")}


## CV LOOP
############

#Some initializations
LambdaRange <- sort(LambdaRange)
hRange <- sort(hRange)
ntrainCV <- dim(Xtrain)[1]-1
nL <- length(LambdaRange)
nh <- length(hRange)
ResCV <- matrix(0,nrow=nL,ncol=nh)
CVG <- matrix(1,nrow=nL,ncol=nh)

for (ncv in 1:ntrain) {

# Determine the data matrix

cvXtrain <- Xtrain[-ncv,]
cvXtest <- matrix(Xtrain[ncv,],nrow=1)
p <- dim(cvXtrain)[2]
r <- min(p,ntrainCV)
#  Standardize the cvXtrain matrix

Sigma2train <- apply(cvXtrain,2,var)*(ntrainCV-1)/ntrainCV
if (sum(Sigma2train==0)!=0){
    if (sum(Sigma2train==0)>(p-2)){
        stop("Message from mgsim.cv.R: the procedure stops because, after leaving one sample, number of predictor variables with no null variance is less than 1.")}
    cvXtrain <- cvXtrain[,which(Sigma2train!=0)]
    cvXtest <- matrix(cvXtest[,which(Sigma2train!=0)],nrow=1)
    Sigma2train <-Sigma2train[which(Sigma2train!=0)]
    p <- dim(cvXtrain)[2]
    r <- min(p,ntrainCV)}
MeancvXtrain <- apply(cvXtrain,2,mean)
sXtrain <- sweep(cvXtrain,2,MeancvXtrain,FUN="-")
sXtrain <- sweep(sXtrain,2,sqrt(Sigma2train),FUN="/")

# Move in the reduced space when necessary
if (p>ntrainCV)
{svd.sXtrain <- svd(t(sXtrain))
 r<-length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
 V <- svd.sXtrain$u[,1:r]
 D <- diag(c(svd.sXtrain$d[1:r]))
 U <- svd.sXtrain$v[,1:r]
 sXtrain <- U%*%D
 rm(D)
 rm(U)
 rm(svd.sXtrain)}
sXtest <- sweep(cvXtest,2,MeancvXtrain,FUN="-")
sXtest <- sweep(sXtest,2,sqrt(Sigma2train),FUN="/")
if (p>ntrainCV)
{sXtest <- sXtest%*%V
 rm(V)
}
rm(cvXtrain)

#Compute Zblock
Z <- cbind(rep(1,ntrainCV),sXtrain)
Zt <- cbind(rep(1,1),sXtest)
Zbloc <- matrix(0,nrow=ntrainCV*c,ncol=c*(r+1))
Ztestbloc <- matrix(0,nrow=c,ncol=c*(r+1))
for (cc in 1:c) {
  row <- (0:(ntrainCV-1))*c+cc
  col <- (r+1)*(cc-1)+1:(r+1)
  Zbloc[row,col] <- Z
  row <- cc
  Ztestbloc[row,col] <- Zt
}
rm(Z)
rm(Zt)


for (j in 1:nh) {
# Compute the Kernel matrix
 WKernel <- matrix(0,ntrainCV,ntrainCV)
 for (jj in 1:ntrainCV) {
    WKernel[,jj] <- apply(exp(-(sXtrain-rep(1,ntrainCV)%*%t(sXtrain[jj,]))^2/(2*hRange[j]^2)),1,prod)
    WKernel[,jj] <- WKernel[,jj]/sum(WKernel[,jj])
    }

 for (i in 1:nL) {
   if (CVG[i,j]==1) {
    res <- mgsimaux(Ytrain[-ncv],Zbloc,sXtrain,LambdaRange[i],WKernel,Ztestbloc,NbIterMax=NbIterMax)
    if (res$Convergence==0)
     {CVG[i,j] <- 0
      ResCV[i,j] <- ntrain}
    if (res$Convergence==1)
     {ResCV[i,j] <- ResCV[i,j]+abs(res$hatY-Ytrain[ncv])}
   }
  }
 }
}


## CONCLUDE
##############

if (sum(CVG)==0)
{stop("No optimal Lambda and h for the given Range")}

#else
# Determine optimal Lambda and h
aux <- which.min(ResCV)
if (nh==1) {
  Lambda <- LambdaRange[aux]
  h <- hRange[1]}

if (nh!=1) {
  h <- hRange[(aux-1)%/%nL+1]
  Lambda <- LambdaRange[(aux-1)%%nL+1]}

return(list(Lambda=Lambda,h=h))

}

