rsurGibbs=
function(Data,Prior,Mcmc)
{
# 
# revision history:
#          P. Rossi 9/05
#          3/07 added classes
#          9/14 changed to improve computations by avoiding Kronecker products
#          W. Taylor 4/15 - added nprint option to MCMC argument
# Purpose:
#   implement Gibbs Sampler for SUR
# 
# Arguments:
#   Data -- regdata
#           regdata is a list of lists of data for each regression
#           regdata[[i]] contains data for regression equation i
#           regdata[[i]]$y is y, regdata[[i]]$X is X
#           note: each regression can have differing numbers of X vars
#                 but you must have same no of obs in each equation. 
#   Prior -- list of prior hyperparameters
#     betabar,A      prior mean, prior precision
#     nu, V          prior on Sigma
#   Mcmc -- list of MCMC parms
#     R number of draws
#     keep -- thinning parameter
#     nprint - print estimated time remaining on every nprint'th draw
# 
# Output: 
#   list of betadraw,Sigmadraw
#
# Model:
#   y_i = X_ibeta + e_i  
#          y is nobs x 1
#          X is nobs x k_i
#          beta is k_i x 1 vector of coefficients
#          i=1,nreg total regressions
#
#         (e_1,k,...,e_nreg,k) ~ N(0,Sigma) k=1,...,nobs
#
#   we can also write as stacked regression
#   y = Xbeta+e
#       y is nobs*nreg x 1,X is nobs*nreg x (sum(k_i))
#   routine draws beta -- the stacked vector of all coefficients
#
# Priors:  beta ~ N(betabar,A^-1)
#          Sigma ~ IW(nu,V)
# 
#
# check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata")}
    if(is.null(Data$regdata)) {pandterm("Requires Data element regdata")}
    regdata=Data$regdata
#
# check regdata for validity
#
nreg=length(regdata)
nobs=length(regdata[[1]]$y)
nvar=0
indreg=double(nreg+1)
y=NULL
for (reg in 1:nreg) {
   if(length(regdata[[reg]]$y) != nobs || nrow(regdata[[reg]]$X) != nobs)
      {pandterm(paste("incorrect dimensions for regression",reg))}
   else
      {indreg[reg]=nvar+1
       nvar=nvar+ncol(regdata[[reg]]$X); y=c(y,regdata[[reg]]$y)}
} 
indreg[nreg+1]=nvar+1
#
# check for Prior
#
if(missing(Prior))
   { betabar=c(rep(0,nvar)); A=BayesmConstant.A*diag(nvar); nu=nreg+BayesmConstant.nuInc; V=nu*diag(nreg)}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))} 
       else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=BayesmConstant.A*diag(nvar)} 
       else {A=Prior$A}
    if(is.null(Prior$nu)) {nu=nreg+BayesmConstant.nuInc}
       else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(nreg)}
       else {ssq=Prior$V}
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
    if(is.null(Mcmc$R)) {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
    if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
    if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
      if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
   }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Gibbs Sampler for SUR Regression Model",fill=TRUE)
cat("  with ",nreg," regressions",fill=TRUE)
cat("  and  ",nobs," observations for each regression",fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat("nu = ",nu,fill=TRUE)
cat("V = ",fill=TRUE)
print(V)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat(" ",fill=TRUE)

#
# set initial value of Sigma
#
E=matrix(double(nobs*nreg),ncol=nreg)
for (reg in 1:nreg) {
    E[,reg]=lm(y~.-1,data=data.frame(y=regdata[[reg]]$y,regdata[[reg]]$X))$residuals
}
Sigma=(crossprod(E)+diag(.01,nreg))/nobs
Sigmainv=chol2inv(chol(Sigma))

#
# precompute various moments and indices into moment matrix and Abetabar
nk=integer(nreg)
Xstar=NULL
Y=NULL
for(i in 1:nreg){
  nk[i]=ncol(regdata[[i]]$X)
  Xstar=cbind(Xstar,regdata[[i]]$X)
  Y=cbind(Y,regdata[[i]]$y)
}
cumnk=cumsum(nk)
XspXs=crossprod(Xstar)
Abetabar=A%*%betabar

###################################################################
# Keunwoo Kim
# 09/19/2014
###################################################################
draws=rsurGibbs_rcpp_loop(regdata,indreg,cumnk,nk,XspXs,Sigmainv,A,Abetabar,nu,V,nvar,E,Y,R,keep,nprint)
###################################################################

attributes(draws$betadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$betadraw)$mcpar=c(1,R,keep)
attributes(draws$Sigmadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(draws$Sigmadraw)$mcpar=c(1,R,keep)

return(draws)
}
