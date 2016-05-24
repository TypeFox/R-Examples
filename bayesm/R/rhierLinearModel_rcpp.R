rhierLinearModel=
function(Data,Prior,Mcmc)
{
#
# Revision History
#     1/17/05  P. Rossi
#     10/05  fixed error in setting prior if Prior argument is missing 
#     3/07 added classes
#     W. Taylor 4/15 - added nprint option to MCMC argument
#
# Purpose:
#   run hiearchical regression model
#
# Arguments:
#   Data list of regdata,Z 
#     regdata is a list of lists each list with members y, X
#        e.g. regdata[[i]]=list(y=y,X=X)
#     X has nvar columns
#     Z is nreg=length(regdata) x nz
#   Prior list of prior hyperparameters
#     Deltabar,A, nu.e,ssq,nu,V
#          note: ssq is a nreg x 1 vector!
#   Mcmc
#     list of Mcmc parameters
#     R is number of draws
#     keep is thining parameter -- keep every keepth draw
#     nprint - print estimated time remaining on every nprint'th draw
#
# Output: 
#   list of 
#   betadraw -- nreg x nvar x R/keep array of individual regression betas
#   taudraw -- R/keep x nreg  array of error variances for each regression
#   Deltadraw -- R/keep x nz x nvar array of Delta draws
#   Vbetadraw -- R/keep x nvar*nvar array of Vbeta draws
#
# Model:
# nreg regression equations 
#        y_i = X_ibeta_i + epsilon_i  
#        epsilon_i ~ N(0,tau_i)
#             nvar X vars in each equation
#
# Priors:
#        tau_i ~ nu.e*ssq_i/chisq(nu.e)  tau_i is the variance of epsilon_i
#        beta_i ~ N(ZDelta[i,],V_beta)
#               Note:  ZDelta is the matrix Z * Delta; [i,] refers to ith row of this product!
#
#          vec(Delta) | V_beta ~ N(vec(Deltabar),Vbeta (x) A^-1)
#          V_beta ~ IW(nu,V)  or V_beta^-1 ~ W(nu,V^-1)
#              Delta, Deltabar are nz x nvar
#              A is nz x nz
#              Vbeta is nvar x nvar
#        
#          NOTE: if you don't have any z vars, set Z=iota (nreg x 1) 
#
#
#  create needed functions
#
#------------------------------------------------------------------------------
append=function(l) { l=c(l,list(XpX=crossprod(l$X),Xpy=crossprod(l$X,l$y)))}
#
getvar=function(l) { 
     v=var(l$y)
     if(is.na(v)) return(1)
     if(v>0) return (v) else return (1)}
#
#------------------------------------------------------------------------------
#

#
# check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata and Z")}
    if(is.null(Data$regdata)) {pandterm("Requires Data element regdata")}
    regdata=Data$regdata
nreg=length(regdata)
if(is.null(Data$Z)) { cat("Z not specified -- putting in iota",fill=TRUE); fsh() ; Z=matrix(rep(1,nreg),ncol=1)}
  else {if (nrow(Data$Z) != nreg) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number regressions ",nreg))}
      else {Z=Data$Z}}
nz=ncol(Z)
#
# check data for validity
#
dimfun=function(l) {c(length(l$y),dim(l$X))}
dims=sapply(regdata,dimfun)
dims=t(dims)
nvar=quantile(dims[,3],prob=.5)

for (i in 1:nreg) 
{
   if(dims[i,1] != dims[i,2]  || dims[i,3] !=nvar) 
      {pandterm(paste("Bad Data dimensions for unit ",i," dims(y,X) =",dims[i,]))}
}
#
# check for Prior
#
if(missing(Prior))
   { Deltabar=matrix(rep(0,nz*nvar),ncol=nvar); A=BayesmConstant.A*diag(nz);
     nu.e=BayesmConstant.nu.e; ssq=sapply(regdata,getvar) ; nu=nvar+BayesmConstant.nuInc ; V= nu*diag(nvar)}
else
   {
    if(is.null(Prior$Deltabar)) {Deltabar=matrix(rep(0,nz*nvar),ncol=nvar)} 
       else {Deltabar=Prior$Deltabar}
    if(is.null(Prior$A)) {A=BayesmConstant.A*diag(nz)} 
       else {A=Prior$A}
    if(is.null(Prior$nu.e)) {nu.e=BayesmConstant.nu.e} 
       else {nu.e=Prior$nu.e}
    if(is.null(Prior$ssq)) {ssq=sapply(regdata,getvar)} 
       else {ssq=Prior$ssq}
    if(is.null(Prior$nu)) {nu=nvar+BayesmConstant.nuInc} 
       else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(nvar)} 
       else {V=Prior$V}
   }
#
# check dimensions of Priors
#
if(ncol(A) != nrow(A) || ncol(A) != nz || nrow(A) != nz) 
   {pandterm(paste("bad dimensions for A",dim(A)))}
if(nrow(Deltabar) != nz || ncol(Deltabar) != nvar)
   {pandterm(paste("bad dimensions for Deltabar ",dim(Deltabar)))}
if(length(ssq) != nreg) {pandterm(paste("bad length for ssq ",length(ssq)))}
if(ncol(V) != nvar || nrow(V) != nvar) {pandterm(paste("bad dimensions for V ",dim(V)))}
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
cat("Starting Gibbs Sampler for Linear Hierarchical Model",fill=TRUE)
cat("   ",nreg," Regressions",fill=TRUE)
cat("   ",ncol(Z)," Variables in Z (if 1, then only intercept)",fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("Deltabar",fill=TRUE)
print(Deltabar)
cat("A",fill=TRUE)
print(A)
cat("nu.e (d.f. parm for regression error variances)= ",nu.e,fill=TRUE)
cat("Vbeta ~ IW(nu,V)",fill=TRUE)
cat("nu = ",nu,fill=TRUE)
cat("V ",fill=TRUE)
print(V)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat(" ",fill=TRUE)
#
#  allocate space for the draws and set initial values of Vbeta and Delta
#
tau=double(nreg)
Delta=matrix(0,nz,nvar)
Vbeta=diag(nvar)
#
#  set up fixed parms for the draw of Vbeta,Delta
#
#  note: in the notation of the MVR  Y =    X      B  
#                                  n x m  n x k  k x m
#                           "n" = nreg
#                           "m" = nvar
#                           "k" = nz
#			general model: Beta = Z Delta + U
#
#       Create XpX elements of regdata and initialize tau
#
regdata=lapply(regdata,append)
tau=sapply(regdata,getvar)

###################################################################
# Keunwoo Kim
# 08/20/2014
###################################################################
draws=rhierLinearModel_rcpp_loop(regdata,Z,Deltabar,A,nu,V,nu.e,ssq,tau,Delta,Vbeta,R,keep,nprint)
###################################################################

attributes(draws$taudraw)$class=c("bayesm.mat","mcmc")
attributes(draws$taudraw)$mcpar=c(1,R,keep)
attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
attributes(draws$Deltadraw)$mcpar=c(1,R,keep)
attributes(draws$Vbetadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(draws$Vbetadraw)$mcpar=c(1,R,keep)
attributes(draws$betadraw)$class=c("bayesm.hcoef")

return(draws)
}