rmvpGibbs=function(Data,Prior,Mcmc){
#
# Revision History:
#   modified by rossi 12/18/04 to include error checking
#   3/07 added classes
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
#
# purpose:  Gibbs MVP model with full covariance matrix
#
# Arguments:
#   Data contains 
#      p the number of alternatives (could be time or could be from pick j of p survey)
#      y -- a vector of length n*p of indicators (1 if "chosen" if not)
#      X -- np x k matrix of covariates (including intercepts)
#                 each X_i is p x nvar
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
#    w_i = X_ibeta + e    e~N(0,Sigma)     note w_i,e are p x 1
#    y_ij = 1 if w_ij > 0 else y_ij = 0  
#  
#  priors:
#    beta ~ N(betabar,A^-1) in prior
#    Sigma ~ IW(nu,V)
#
#  Check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of p, y, X")}
if(is.null(Data$p)) {pandterm("Requires Data element p -- number of binary indicators")}
p=Data$p
if(is.null(Data$y)) {pandterm("Requires Data element y -- values of binary indicators")}
y=Data$y
if(is.null(Data$X)) {pandterm("Requires Data element X -- matrix of covariates")}
X=Data$X
#
# check data for validity
#
levely=as.numeric(levels(as.factor(y)))
bady=FALSE
for (i in 0:1) 
{ if(levely[i+1] != i) {bady=TRUE} }
cat("Table of y values",fill=TRUE)
print(table(y))
if (bady) {pandterm("Invalid y")}
if (length(y)%%p !=0) {pandterm("length of y is not a multiple of p")}
n=length(y)/p
k=ncol(X)
if(nrow(X) != (n*p)) {pandterm(paste("X has ",nrow(X)," rows; must be = p*n"))}
#
# check for prior elements
#
if(missing(Prior)) 
{ betabar=rep(0,k) ; A=BayesmConstant.A*diag(k) ; nu=p+BayesmConstant.nuInc; V=nu*diag(p)}
else 
{if(is.null(Prior$betabar)) {betabar=rep(0,k)} else {betabar=Prior$betabar}
 if(is.null(Prior$A)) {A=BayesmConstant.A*diag(k)} else {A=Prior$A}
 if(is.null(Prior$nu)) {nu=p+BayesmConstant.nuInc} else {nu=Prior$nu}
 if(is.null(Prior$V)) {V=nu*diag(p)} else {V=Prior$V}}
if(length(betabar) != k) pandterm("length betabar ne k")
if(sum(dim(A)==c(k,k)) != 2) pandterm("A is of incorrect dimension")
if(nu < 1) pandterm("invalid nu value")
if(sum(dim(V)==c(p,p)) != 2) pandterm("V is of incorrect dimension")
#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R must be included")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$beta0)) {beta0=rep(0,k)} else {beta0=Mcmc$beta0}
if(is.null(Mcmc$sigma0)) {sigma0=diag(p)} else {sigma0=Mcmc$sigma0}
if(length(beta0) != k) pandterm("beta0 is not of length k")
if(sum(dim(sigma0) == c(p,p)) != 2) pandterm("sigma0 is of incorrect dimension")
if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
  if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting Gibbs Sampler for MVP",fill=TRUE)
cat("  ",n," obs of ",p," binary indicators; ",k," indep vars (including intercepts)",fill=TRUE)
cat(" ",fill=TRUE)
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
cat("initial sigma= ",fill=TRUE)
print(sigma0)
cat(" ",fill=TRUE)

###################################################################
# Wayne Taylor
# 09/03/2014
###################################################################
loopout = rmvpGibbs_rcpp_loop(R,keep,nprint,p,y,X,beta0,sigma0,V,nu,betabar,A);
###################################################################

attributes(loopout$betadraw)$class=c("bayesm.mat","mcmc")
attributes(loopout$betadraw)$mcpar=c(1,R,keep)
attributes(loopout$sigmadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(loopout$sigmadraw)$mcpar=c(1,R,keep)

return(loopout)
}