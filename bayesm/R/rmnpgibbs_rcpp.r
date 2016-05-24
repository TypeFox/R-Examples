rmnpGibbs=function(Data,Prior,Mcmc) {
#
# Revision History:
#   modified by rossi 12/18/04 to include error checking
#   3/07 added classes
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
# purpose:  Gibbs MNP model with full covariance matrix
#
# Arguments:
#   Data contains 
#      p the number of choice alternatives
#      y -- a vector of length n with choices (takes on values from 1, .., p)
#      X -- n(p-1) x k matrix of covariates (including intercepts)
#           note: X is the differenced matrix unlike MNL X=stack(X_1,..,X_n) 
#                 each X_i is (p-1) x nvar
#
#   Prior contains a list of (betabar, A, nu, V)
#      if elements of prior do not exist, defaults are used
#
#   Mcmc is a list of (beta0,sigma0,R,keep)  
#     beta0,sigma0 are intial values, if not supplied defaults are used
#     R is number of draws
#     keep is thinning parm, keep every keepth draw
#     nprint - print estimated time remaining on every nprint'th draw
#
# Output: a list of every keepth betadraw and sigmsdraw
#
#  model: 
#    w_i = X_ibeta + e    e~N(0,Sigma)     note w_i,e are (p-1) x 1
#    y_i = j  if w_ij > w_i-j  j=1,...,p-1
#    y_i = p  if all w_i < 0
#  
#  priors:
#    beta ~ N(betabar,A^-1)
#    Sigma ~ IW(nu,V)
#
#  Check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of p, y, X")}
  if(is.null(Data$p)) {pandterm("Requires Data element p -- number of alternatives")}
  p=Data$p
  if(is.null(Data$y)) {pandterm("Requires Data element y -- number of alternatives")}
  y=Data$y
  if(is.null(Data$X)) {pandterm("Requires Data element X -- matrix of covariates")}
  X=Data$X
#
# check data for validity
#
levely=as.numeric(levels(as.factor(y)))
if(length(levely) != p) {pandterm(paste("y takes on ",length(levely),
  " values -- must be ",p))}
  bady=FALSE
  for (i in 1:p) 
  {
      if(levely[i] != i) bady=TRUE
  }
cat("Table of y values",fill=TRUE)
print(table(y))
if (bady) {pandterm("Invalid y")}
n=length(y)
k=ncol(X)
pm1=p-1
if(nrow(X)/n != pm1) {pandterm(paste("X has ",nrow(X)," rows; must be = (p-1)n"))}
#
# check for prior elements
#
if(missing(Prior)) 
  { betabar=rep(0,k) ; A=BayesmConstant.A*diag(k) ; nu=pm1+3; V=nu*diag(pm1)}
else 
  {if(is.null(Prior$betabar)) {betabar=rep(0,k)} else {betabar=Prior$betabar}
   if(is.null(Prior$A)) {A=BayesmConstant.A*diag(k)} else {A=Prior$A}
   if(is.null(Prior$nu)) {nu=pm1+BayesmConstant.nuInc} else {nu=Prior$nu}
   if(is.null(Prior$V)) {V=nu*diag(pm1)} else {V=Prior$V}}
if(length(betabar) != k) pandterm("length betabar ne k")
if(sum(dim(A)==c(k,k)) != 2) pandterm("A is of incorrect dimension")
if(nu < 1) pandterm("invalid nu value")
if(sum(dim(V)==c(pm1,pm1)) != 2) pandterm("V is of incorrect dimension")
#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R must be included")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$beta0)) {beta0=rep(0,k)} else {beta0=Mcmc$beta0}
if(is.null(Mcmc$sigma0)) {sigma0=diag(pm1)} else {sigma0=Mcmc$sigma0}
if(length(beta0) != k) pandterm("beta0 is not of length k")
if(sum(dim(sigma0) == c(pm1,pm1)) != 2) pandterm("sigma0 is of incorrect dimension")
if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
  if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting Gibbs Sampler for MNP",fill=TRUE)
cat("  ",n," obs; ",p," choice alternatives; ",k," indep vars (including intercepts)",fill=TRUE)
cat(" ",fill=TRUE)
cat("Table of y values",fill=TRUE)
print(table(y))
cat("Prior Parms:",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat("nu",fill=TRUE)
print(nu)
cat("V",fill=TRUE)
print(V)
cat(" ",fill=TRUE)
cat("MCMC Parms:",fill=TRUE)
cat("  ",R," reps; keeping every ",keep,"th draw"," nprint= ",nprint,fill=TRUE)
cat("initial beta= ",beta0,fill=TRUE)
cat("initial sigma= ",sigma0,fill=TRUE)
cat(" ",fill=TRUE)

###################################################################
# Wayne Taylor
# 09/03/2014
###################################################################
loopout = rmnpGibbs_rcpp_loop(R,keep,nprint,pm1,y,X,beta0,sigma0,V,nu,betabar,A);
###################################################################

attributes(loopout$betadraw)$class=c("bayesm.mat","mcmc")
attributes(loopout$betadraw)$mcpar=c(1,R,keep)
attributes(loopout$sigmadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(loopout$sigmadraw)$mcpar=c(1,R,keep)

return(loopout)
}
