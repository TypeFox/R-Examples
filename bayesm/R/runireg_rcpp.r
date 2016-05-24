runireg=
function(Data,Prior,Mcmc)
{
# 
# revision history:
#          P. Rossi 1/17/05
#          revised 9/05 to put in Data,Prior,Mcmc calling convention
#          3/07 added classes
#          W. Taylor 4/15 - added nprint option to MCMC argument
# Purpose:
#   perform iid draws from posterior of regression model using
#     conjugate prior
# 
# Arguments:
#   Data -- list of data 
#           y,X
#   Prior -- list of prior hyperparameters
#     betabar,A      prior mean, prior precision
#     nu, ssq        prior on sigmasq
#   Mcmc -- list of MCMC parms
#     R number of draws
#     keep -- thinning parameter
#     nprint - print estimated time remaining on every nprint'th draw
# 
# Output: 
#   list of beta, sigmasq
#
# Model:
#   y = Xbeta + e  e ~N(0,sigmasq)
#          y is n x 1
#          X is n x k
#          beta is k x 1 vector of coefficients
#
# Priors:  beta ~ N(betabar,sigmasq*A^-1)
#          sigmasq ~ (nu*ssq)/chisq_nu
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
if(nobs != nrow(X) ) {pandterm("length(y) ne nrow(X)")}
#
# check for Prior
#
if(missing(Prior))
   { betabar=c(rep(0,nvar)); A=BayesmConstant.A*diag(nvar); nu=BayesmConstant.nu; ssq=var(y)}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))} 
       else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=BayesmConstant.A*diag(nvar)} 
       else {A=Prior$A}
    if(is.null(Prior$nu)) {nu=BayesmConstant.nu}
       else {nu=Prior$nu}
    if(is.null(Prior$ssq)) {ssq=var(y)}
       else {ssq=Prior$ssq}
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
cat("Starting IID Sampler for Univariate Regression Model",fill=TRUE)
cat("  with ",nobs," observations",fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat("nu = ",nu," ssq= ",ssq,fill=TRUE)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat(" ",fill=TRUE)

###################################################################
# Keunwoo Kim
# 08/05/2014
###################################################################
draws = runireg_rcpp_loop(y, X, betabar, A, nu, ssq, R, keep, nprint)
###################################################################

attributes(draws$betadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$betadraw)$mcpar=c(1,R,keep)
attributes(draws$sigmasqdraw)$class=c("bayesm.mat","mcmc")
attributes(draws$sigmasqdraw)$mcpar=c(1,R,keep)

return(draws)
}
