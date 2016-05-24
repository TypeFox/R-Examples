rmnlIndepMetrop=function(Data,Prior,Mcmc){
#
# revision history:
#   p. rossi 1/05
#   2/9/05 fixed error in Metrop eval
#   changed to reflect new argument order in llmnl,mnlHess 9/05
#   added return for log-like  11/05
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
# purpose: 
#   draw from posterior for MNL using Independence Metropolis
#
# Arguments:
#   Data - list of p,y,X  
#     p is number of alternatives
#     X is nobs*p x nvar matrix
#     y is nobs vector of values from 1 to p
#   Prior - list of A, betabar
#     A is nvar x nvar prior preci matrix
#     betabar is nvar x 1 prior mean
#   Mcmc
#     R is number of draws
#     keep is thinning parameter
#     nprint - print estimated time remaining on every nprint'th draw
#     nu degrees of freedom parameter for independence 
#        sampling density
#
# Output:
#   list of betadraws
#
# Model:   Pr(y=j) = exp(x_j'beta)/sum(exp(x_k'beta)
#
# Prior:   beta ~ N(betabar,A^-1)
#
# check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of p, y, X")}
    if(is.null(Data$X)) {pandterm("Requires Data element X")}
    X=Data$X
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=Data$y
    if(is.null(Data$p)) {pandterm("Requires Data element p")}
    p=Data$p
nvar=ncol(X)
nobs=length(y)
#
# check data for validity
#
if(length(y) != (nrow(X)/p) ) {pandterm("length(y) ne nrow(X)/p")}
if(sum(y %in% (1:p)) < nobs) {pandterm("invalid values in y vector -- must be integers in 1:p")}
cat(" table of y values",fill=TRUE)
print(table(y))
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
    if(is.null(Mcmc$nu)) {nu=6} else {nu=Mcmc$nu}
   }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Independence Metropolis Sampler for Multinomial Logit Model",fill=TRUE)
cat("  ",length(y)," obs with ",p," alternatives",fill=TRUE)
cat(" ", fill=TRUE)
cat("Table of y Values",fill=TRUE)
print(table(y))
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint," nu (df for st candidates) = ",nu,fill=TRUE)
cat(" ",fill=TRUE)

#
# compute required quantities for indep candidates
#
beta=c(rep(0,nvar))
mle=optim(beta,llmnl,X=X,y=y,method="BFGS",hessian=TRUE,control=list(fnscale=-1))
beta=mle$par
betastar=mle$par
mhess=mnlHess(beta,y,X)
candcov=chol2inv(chol(mhess))
root=chol(candcov)
rooti=backsolve(root,diag(nvar))
priorcov=chol2inv(chol(A))
rootp=chol(priorcov)
rootpi=backsolve(rootp,diag(nvar))

oldloglike=llmnl(beta,y,X)
oldlpost=oldloglike+lndMvn(beta,betabar,rootpi)
oldlimp=lndMvst(beta,nu,betastar,rooti)
#       note: we don't need the determinants as they cancel in
#       computation of acceptance prob

###################################################################
# Wayne Taylor
# 08/21/2014
###################################################################
loopout = rmnlIndepMetrop_rcpp_loop(R,keep,nu,betastar,root,y,X,betabar,rootpi,rooti,oldlimp,oldlpost,nprint);
###################################################################

attributes(loopout$betadraw)$class=c("bayesm.mat","mcmc")
attributes(loopout$betadraw)$mcpar=c(1,R,keep)

return(list(betadraw=loopout$betadraw,loglike=loopout$loglike,acceptr=loopout$naccept/R))
}
