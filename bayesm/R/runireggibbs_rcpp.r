runiregGibbs=
function(Data,Prior,Mcmc)
{
# 
# revision history:
#          P. Rossi 1/17/05
#          3/07 added classes
#          W. Taylor 4/15 - added nprint option to MCMC argument
# Purpose:
#   perform Gibbs iterations for Univ Regression Model using
#     prior with beta, sigma-sq indep
# 
# Arguments:
#   Data -- list of data 
#           y,X
#   Prior -- list of prior hyperparameters
#     betabar,A      prior mean, prior precision
#     nu, ssq        prior on sigmasq
#   Mcmc -- list of MCMC parms
#     sigmasq=initial value for sigmasq
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
# Priors:  beta ~ N(betabar,A^-1)
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
   { betabar=c(rep(0,nvar)); A=.01*diag(nvar); nu=3; ssq=var(y)}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))} 
       else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=.01*diag(nvar)} 
       else {A=Prior$A}
    if(is.null(Prior$nu)) {nu=3}
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
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    if(is.null(Mcmc$nprint)) {nprint=100} else {nprint=Mcmc$nprint}
      if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
    if(is.null(Mcmc$sigmasq)) {sigmasq=var(y)} else {sigmasq=Mcmc$sigmasq}
   }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Gibbs Sampler for Univariate Regression Model",fill=TRUE)
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
draws = runiregGibbs_rcpp_loop(y, X, betabar, A, nu, ssq, sigmasq, R, keep, nprint)
###################################################################

attributes(draws$betadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$betadraw)$mcpar=c(1,R,keep)
attributes(draws$sigmasqdraw)$class=c("bayesm.mat","mcmc")
attributes(draws$sigmasqdraw)$mcpar=c(1,R,keep)

return(draws)
}
