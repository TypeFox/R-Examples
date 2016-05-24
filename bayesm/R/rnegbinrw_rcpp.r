rnegbinRw=function(Data, Prior, Mcmc){
#   Revision History
#	  Sridhar Narayanan - 05/2005
#         P. Rossi 6/05
#         3/07 added classes
#   Keunwoo Kim 11/2014
#         1. added "alphafix" argument 
#         2. corrected code in more clear way (in Cpp)
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
#   Model
#       (y|lambda,alpha) ~ Negative Binomial(Mean = lambda, Overdispersion par = alpha)
#
#       ln(lambda) =  X * beta
#               
#   Priors
#       beta ~ N(betabar, A^-1)
#       alpha ~ Gamma(a,b) where mean = a/b and variance = a/(b^2)
#
#   Arguments
#       Data = list of y, X
#              e.g. regdata[[i]]=list(y=y,X=X)
#              X has nvar columns including a first column of ones
#
#       Prior - list containing the prior parameters
#           betabar, A - mean of beta prior, inverse of variance covariance of beta prior
#           a, b - parameters of alpha prior
#
#       Mcmc - list containing
#           R is number of draws
#           keep is thinning parameter (def = 1)
#           nprint - print estimated time remaining on every nprint'th draw (def = 100)
#           s_beta - scaling parameter for beta RW (def = 2.93/sqrt(nvar))
#           s_alpha - scaling parameter for alpha RW (def = 2.93)
#           beta0 - initial guesses for parameters, if not supplied default values are used
#           alpha - value of alpha fixed. If it is given, draw only beta
#

#
# Definitions of functions used within rhierNegbinRw
#
llnegbin = 
function(par,X,y, nvar) {
# Computes the log-likelihood
    beta = par[1:nvar]
    alpha = exp(par[nvar+1])+1.0e-50
    mean=exp(X%*%beta)
    prob=alpha/(alpha+mean)
    prob=ifelse(prob<1.0e-100,1.0e-100,prob)
     out=dnbinom(y,size=alpha,prob=prob,log=TRUE)
     return(sum(out))
}

#
# Error Checking
#
if(missing(Data)) {pandterm("Requires Data argument -- list of X and y")}
if(is.null(Data$X)) {pandterm("Requires Data element X")} else {X=Data$X}
if(is.null(Data$y)) {pandterm("Requires Data element y")} else {y=Data$y}
nvar = ncol(X)

if (length(y) != nrow(X)) {pandterm("Mismatch in the number of observations in X and y")}
nobs=length(y)

#
# check for prior elements
#
if(missing(Prior)) {
    betabar=rep(0,nvar); A=BayesmConstant.A*diag(nvar) ;  a=BayesmConstant.agammaprior; b=BayesmConstant.bgammaprior;
}
else {
    if(is.null(Prior$betabar)) {betabar=rep(0,nvar)} else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=BayesmConstant.A*diag(nvar)} else {A=Prior$A}
    if(is.null(Prior$a)) {a=BayesmConstant.agammaprior} else {a=Prior$a}
    if(is.null(Prior$b)) {b=BayesmConstant.bgammaprior} else {b=Prior$b}
}

if(length(betabar) != nvar) pandterm("betabar is of incorrect dimension")
if(sum(dim(A)==c(nvar,nvar)) != 2) pandterm("A is of incorrect dimension")
if((length(a) != 1) | (a <=0)) pandterm("a should be a positive number")
if((length(b) != 1) | (b <=0)) pandterm("b should be a positive number")

#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$beta0)) {beta0=rep(0,nvar)} else {beta0=Mcmc$beta0}
if(length(beta0) !=nvar) pandterm("beta0 is not of dimension nvar")
if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
  if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
if(is.null(Mcmc$s_alpha)) {cat("Using default s_alpha = 2.93",fill=TRUE); s_alpha=BayesmConstant.RRScaling} 
    else {s_alpha = Mcmc$s_alpha} 
if(is.null(Mcmc$s_beta)) {cat("Using default s_beta = 2.93/sqrt(nvar)",fill=TRUE); s_beta=BayesmConstant.RRScaling/sqrt(nvar)} 
    else {s_beta = Mcmc$s_beta}
# Keunwoo Kim 11/2014 #############################################
if(is.null(Mcmc$alpha)) {fixalpha=FALSE} else {fixalpha=TRUE; alpha=Mcmc$alpha}
if(fixalpha & alpha<=0) pandterm("alpha is not positive")
###################################################################

#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting Random Walk Metropolis Sampler for Negative Binomial Regression",fill=TRUE)
cat("  ",nobs," obs; ",nvar," covariates (including intercept); ",fill=TRUE)
cat("Prior Parameters:",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat("a",fill=TRUE)
print(a)
cat("b",fill=TRUE)
print(b)
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat("s_alpha = ",s_alpha,fill=TRUE)
cat("s_beta = ",s_beta,fill=TRUE)
cat(" ",fill=TRUE)

par = rep(0,(nvar+1))
cat(" Initializing RW Increment Covariance Matrix...",fill=TRUE)
fsh()
mle = optim(par,llnegbin, X=X, y=y, nvar=nvar, method="L-BFGS-B", upper=c(Inf,Inf,Inf,log(100000000)), hessian=TRUE, control=list(fnscale=-1))
fsh()
beta_mle=mle$par[1:nvar]
alpha_mle = exp(mle$par[nvar+1])
varcovinv = -mle$hessian
beta = beta0
betacvar = s_beta*solve(varcovinv[1:nvar,1:nvar])
betaroot = t(chol(betacvar))
if(!fixalpha) {alpha = alpha_mle}
alphacvar = s_alpha/varcovinv[nvar+1,nvar+1]
alphacroot = sqrt(alphacvar)
cat("beta_mle = ",beta_mle,fill=TRUE)
cat("alpha_mle = ",alpha_mle, fill = TRUE)
fsh()

###################################################################
# Keunwoo Kim
# 09/03/2014
###################################################################
if (fixalpha) {alpha=Mcmc$alpha}
draws=rnegbinRw_rcpp_loop(y, X, betabar, chol(A), a, b, beta, alpha, fixalpha, betaroot, alphacroot, R, keep, nprint)
###################################################################

attributes(draws$betadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$betadraw)$mcpar=c(1,R,keep)
attributes(draws$alphadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$alphadraw)$mcpar=c(1,R,keep)
return(list(betadraw=draws$betadraw,alphadraw=draws$alphadraw,
            acceptrbeta=draws$nacceptbeta/R*keep,acceptralpha=draws$nacceptalpha/R*keep))
}
