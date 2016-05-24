rbprobitGibbs=
function(Data,Prior,Mcmc)
{
#
# revision history:
#   p. rossi 1/05
#   3/07 added validity check of values of y and classes
#   3/07 fixed error with betabar supplied
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
# purpose: 
#   draw from posterior for binary probit using Gibbs Sampler
#
# Arguments:
#   Data - list of X,y  
#     X is nobs x nvar, y is nobs vector of 0,1
#   Prior - list of A, betabar
#     A is nvar x nvar prior preci matrix
#     betabar is nvar x 1 prior mean
#   Mcmc
#     R is number of draws
#     keep is thinning parameter
#     nprint - print estimated time remaining on every nprint'th draw
#
# Output:
#   list of betadraws
#
# Model:   y = 1 if  w=Xbeta + e   > 0  e ~N(0,1)
#
# Prior:   beta ~ N(betabar,A^-1)
#
#
#
# check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of y and X")}
    if(is.null(Data$X)) {pandterm("Requires Data element X")}
    X=Data$X
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=Data$y
nvar=ncol(X)
nobs=length(y)
#
# check data for validity
#
if(length(y) != nrow(X) ) {pandterm("y and X not of same row dim")}
if(sum(unique(y) %in% c(0:1)) < length(unique(y))) {pandterm("Invalid y, must be 0,1")}
#
# check for Prior
#
if(missing(Prior))
   { betabar=c(rep(0,nvar)); A=BayesmConstant.A*diag(nvar)}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))} 
       else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=BayesmConstant.A*diag(nvar)} 
       else {A=Prior$A}
   }
#
# check dimensions of Priors
#
if(ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) 
   {pandterm(paste("bad dimensions for A",dim(A)))}
if(length(betabar) != nvar)
   {pandterm(paste("betabar wrong length, length= ",length(betabar)))}
#
# check MCMC argument
#
if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
else
   {
    if(is.null(Mcmc$R)) 
       {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
    if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
    if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
      if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
    }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Gibbs Sampler for Binary Probit Model",fill=TRUE)
cat("   with ",length(y)," observations",fill=TRUE)
cat("Table of y Values",fill=TRUE)
print(table(y))
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat(" ",fill=TRUE)

beta=c(rep(0,nvar))
sigma=c(rep(1,nrow(X)))
root=chol(chol2inv(chol((crossprod(X,X)+A))))
Abetabar=crossprod(A,betabar)
        a=ifelse(y == 0,-100, 0)
        b=ifelse(y == 0, 0, 100)

###################################################################
# Keunwoo Kim
# 08/05/2014
###################################################################
draws=rbprobitGibbs_rcpp_loop(y,X,Abetabar,root,beta,sigma,a,b,R,keep,nprint)
###################################################################

attributes(draws$betadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$betadraw)$mcpar=c(1,R,keep)
return(draws)
}
