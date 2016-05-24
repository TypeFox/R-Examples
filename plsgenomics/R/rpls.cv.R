### rpls.cv.R  (2006-01)
###
###    Determination by Cross-validation of Ridge Partial Least square
###                 hyper-parameters for binary data
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

rpls.cv <- function (Ytrain,Xtrain,LambdaRange,ncompMax,NbIterMax=50)
{

##    INPUT VARIABLES
#########################
##  Xtrain   : matrix ntrain x p
##      train data matrix
##  Ytrain   : vector ntrain
##      response variable {1,2}-valued vector
##  LambdaRange : vector nLambda
##      possible values for the regularization parameter Lambda
##  NbIterMax : positive integer
##      maximal number of iteration in the WIRRLS part
##  ncompMax : positive integer
##      maximal number of PLS components
##      ncompMax=0 : Ridge


##    OUTPUT VARIABLES
##########################
## Lambda : optimal regularization parameter Lambda
## ncomp : optimal number of PLS components


##  TEST ON INPUT VARIABLES
##############################
#On Xtrain
if ((is.matrix(Xtrain)==FALSE)||(is.numeric(Xtrain)==FALSE)) {
 stop("Message from rpls.cv.R: Xtrain is not of valid type")}

if (dim(Xtrain)[2]==1) {
 stop("Message from rpls.cv.R: p=1 is not valid")}

ntrain <- dim(Xtrain)[1]

#On Ytrain
if ((is.vector(Ytrain)==FALSE)||(is.numeric(Ytrain)==FALSE)) {
 stop("Message from rpls.cv.R: Ytrain is not of valid type")}

if (length(Ytrain)!=ntrain) {
 stop("Message from rpls.cv.R: the length of Ytrain is not equal to the Xtrain row number")}

Ytrain <- Ytrain-1

if ((sum(floor(Ytrain)-Ytrain)!=0)||(sum(Ytrain<0)>0)){
 stop("Message from rpls.cv.R: Ytrain is not of valid type")}

c <- max(Ytrain)
if (c!=1) {
 stop("Message from rpls.cv.R: Ytrain is not of valid type")}

eff<-rep(0,2)
for (i in 0:1) {
    eff[i+1]<-sum(Ytrain==i)}
if (sum(eff<=1)>0) {
 stop("Message from rpls.cv.R: there are not enough samples for each class")}

#On hyper parameters range

if ((is.numeric(LambdaRange)==FALSE)||(is.vector(LambdaRange)==FALSE)||(sum(LambdaRange<0)>0)){
 stop("Message from rpls.cv.R: LambdaRange is not of valid type")}

if ((is.numeric(ncompMax)==FALSE)||(round(ncompMax)-ncompMax!=0)||(ncompMax<0)){
 stop("Message from rpls.cv.R: ncompMax is not of valid type")}

if ((is.numeric(NbIterMax)==FALSE)||(round(NbIterMax)-NbIterMax!=0)||(NbIterMax<1)){
 stop("Message from rpls.cv.R: NbIterMax is not of valid type")}



## CV LOOP
############

#Some initializations
LambdaRange <- sort(LambdaRange)
ntrainCV <- dim(Xtrain)[1]-1


LambdaIndex <- 1:length(LambdaRange)
if (ncompMax <=1){
    nc <- 1}
if (ncompMax > 1){
    nc <- ncompMax}
ResCV <- matrix(0,nrow=length(LambdaRange),ncol=nc)

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
        stop("Message from rpls.cv.R: the procedure stops because, after leaving one sample, number of predictor variables with no null variance is less than 1.")}
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
rm(V)}
rm(cvXtrain)

LIndexaux <- LambdaIndex

for (i in LambdaIndex) {
    res <- rplsaux(Ytrain[-ncv],sXtrain,LambdaRange[i],ncompMax,sXtest,NbIterMax=NbIterMax)
    if (res$Convergence==0)
     {LIndexaux <- LIndexaux[LIndexaux!=i]}
    if (res$Convergence==1)
    {ResCV[i,] <- ResCV[i,]+abs(res$hatY-Ytrain[ncv])}
    }

LambdaIndex <- LIndexaux

}

## CONCLUDE
##############
if (length(LambdaIndex)==0)
{stop("No optimal Lambda for the given LambdaRange")}

#else
ResCV <- ResCV[LambdaIndex,]
if (length(LambdaIndex)==1)  {
    ResCV <- matrix(ResCV,nrow=1)}
# Determine optimal Lambda and ncomp
aux <- which.min(ResCV)
if (ncompMax <=1){
    ncomp <- ncompMax
    Lambda <- LambdaRange[LambdaIndex[aux]]}
if (ncompMax > 1){
    nl <- dim(ResCV)[1]
    ncomp <- (aux-1)%/%nl+1
    Lambda <- LambdaRange[LambdaIndex[(aux-1)%%nl+1]]}

return(list(ncomp=ncomp,Lambda=Lambda))

}

