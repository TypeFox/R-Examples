rhierLinearMixture=function(Data,Prior,Mcmc){
#
# revision history:
#   changed 12/17/04 by rossi to fix bug in drawdelta when there is zero/one unit
#   in a mixture component
#   adapted to linear model by Vicky Chen 6/06
#   put in classes 3/07
#   changed a check 9/08
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
# purpose: run hierarchical linear model with mixture of normals 
#
# Arguments:
#   Data contains a list of (regdata, and possibly Z)
#      regdata is a list of lists (one list per unit)
#          regdata[[i]]=list(y,X)
#             y is a vector of observations
#             X is a length(y) x nvar matrix of values of
#               X vars including intercepts
#             Z is an nreg x nz matrix of values of variables
#               note: Z should NOT contain an intercept
#   Prior contains a list of (nu.e,ssq,deltabar,Ad,mubar,Amu,nu,V,ncomp,a) 
#      ncomp is the number of components in normal mixture
#           if elements of Prior (other than ncomp) do not exist, defaults are used
#   Mcmc contains a list of (s,c,R,keep,nprint)
#
# Output:  as list containing
#   taodraw is R/keep x nreg  array of error variances for each regression
#   Deltadraw R/keep  x nz*nvar matrix of draws of Delta, first row is initial value
#   betadraw is nreg x nvar x R/keep array of draws of betas
#   probdraw is R/keep x ncomp matrix of draws of probs of mixture components
#   compdraw is a list of list of lists (length R/keep)
#      compdraw[[rep]] is the repth draw of components for mixtures
#
# Priors:
#    tau_i ~ nu.e*ssq_i/chisq(nu.e)  tau_i is the variance of epsilon_i
#    beta_i = delta %*% z[i,] + u_i
#       u_i ~ N(mu_ind[i],Sigma_ind[i])
#       ind[i] ~multinomial(p)
#       p ~ dirichlet (a)
#           a: Dirichlet parameters for prior on p
#       delta is a k x nz array
#          delta= vec(D) ~ N(deltabar,A_d^-1)
#    mu_j ~ N(mubar,A_mu^-1(x)Sigma_j)
#    Sigma_j ~ IW(nu,V^-1)
#    ncomp is number of components
#
# MCMC parameters
#   R is number of draws
#   keep is thinning parameter, keep every keepth draw
#   nprint - print estimated time remaining on every nprint'th draw
#
#  check arguments
#
#--------------------------------------------------------------------------------------------------
#
#  create functions needed
#
append=function(l) { l=c(l,list(XpX=crossprod(l$X),Xpy=crossprod(l$X,l$y)))}
#
getvar=function(l) { 
  v=var(l$y)
  if(is.na(v)) return(1)
  if(v>0) return (v) else return (1)}
#
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata, and (possibly) Z")}
  if(is.null(Data$regdata)) {pandterm("Requires Data element regdata (list of data for each unit)")}
  regdata=Data$regdata
  nreg=length(regdata)
  drawdelta=TRUE
if(is.null(Data$Z)) { cat("Z not specified",fill=TRUE); fsh() ; drawdelta=FALSE}
  else {if (nrow(Data$Z) != nreg) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number regressions ",nreg))}
      else {Z=Data$Z}}
  if(drawdelta) {
     nz=ncol(Z)
     colmeans=apply(Z,2,mean)
     if(sum(colmeans) > .00001) 
       {pandterm(paste("Z does not appear to be de-meaned: colmeans= ",colmeans))}
  }
#
# check regdata for validity
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
# check on prior
#
if(missing(Prior)) 
{pandterm("Requires Prior list argument (at least ncomp)")} 
if(is.null(Prior$nu.e)) {nu.e=BayesmConstant.nu.e} 
   else {nu.e=Prior$nu.e}
if(is.null(Prior$ssq)) {ssq=sapply(regdata,getvar)} 
   else {ssq=Prior$ssq}
if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp (num of mixture components)")} else {ncomp=Prior$ncomp}
if(is.null(Prior$mubar)) {mubar=matrix(rep(0,nvar),nrow=1)} else { mubar=matrix(Prior$mubar,nrow=1)}
  if(ncol(mubar) != nvar) {pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ",ncol(mubar)))}
if(is.null(Prior$Amu)) {Amu=matrix(BayesmConstant.A,ncol=1)} else {Amu=matrix(Prior$Amu,ncol=1)}
  if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Am must be a 1 x 1 array")}
if(is.null(Prior$nu)) {nu=nvar+BayesmConstant.nuInc}  else {nu=Prior$nu}
  if(nu < 1) {pandterm("invalid nu value")}
if(is.null(Prior$V)) {V=nu*diag(nvar)} else {V=Prior$V}
  if(sum(dim(V)==c(nvar,nvar)) !=2) pandterm("Invalid V in prior")
if(is.null(Prior$Ad) & drawdelta) {Ad=BayesmConstant.A*diag(nvar*nz)} else {Ad=Prior$Ad}
if(drawdelta) {if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}}
if(is.null(Prior$deltabar)& drawdelta) {deltabar=rep(0,nz*nvar)} else {deltabar=Prior$deltabar}
  if(drawdelta) {if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}}
if(is.null(Prior$a)) { a=rep(BayesmConstant.a,ncomp)} else {a=Prior$a}
if(length(a) != ncomp) {pandterm("Requires dim(a)= ncomp (no of components)")}
bada=FALSE
   for(i in 1:ncomp) { if(a[i] < 0) bada=TRUE}
  if(bada) pandterm("invalid values in a vector")
#
# check on Mcmc
#
if(missing(Mcmc)) 
  {pandterm("Requires Mcmc list argument")}
else 
   { 
    if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {R=Mcmc$R}
    if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
      if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
    }
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting MCMC Inference for Hierarchical Linear Model:",fill=TRUE)
cat("   Normal Mixture with",ncomp,"components for first stage prior",fill=TRUE)
cat(paste("   for ",nreg," cross-sectional units"),fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("nu.e =",nu.e,fill=TRUE)
cat("nu =",nu,fill=TRUE)
cat("V ",fill=TRUE)
print(V)
cat("mubar ",fill=TRUE)
print(mubar)
cat("Amu ", fill=TRUE)
print(Amu)
cat("a ",fill=TRUE)
print(a)
if(drawdelta) 
{
   cat("deltabar",fill=TRUE)
   print(deltabar)
   cat("Ad",fill=TRUE)
   print(Ad)
}
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat("",fill=TRUE)

#  initialize values
#
#  Create XpX elements of regdata and initialize tau
#
regdata=lapply(regdata,append)
tau=sapply(regdata,getvar)
#
# set initial values for the indicators
#     ind is of length(nreg) and indicates which mixture component this obs
#     belongs to.
#
ind=NULL
ninc=floor(nreg/ncomp)
for (i in 1:(ncomp-1)) {ind=c(ind,rep(i,ninc))}
if(ncomp != 1) {ind = c(ind,rep(ncomp,nreg-length(ind)))} else {ind=rep(1,nreg)}
#
#initialize delta
#
if (drawdelta){
  olddelta = rep(0,nz*nvar)
} else { #send placeholders to the _loop function if there is no Z matrix
  olddelta = 0
  Z = matrix(0)
  deltabar = 0
  Ad = matrix(0)
}
#
# initialize probs
#
oldprob=rep(1/ncomp,ncomp)

###################################################################
# Wayne Taylor
# 09/19/2014
###################################################################
draws =  rhierLinearMixture_rcpp_loop(regdata, Z,
                                      deltabar, Ad, mubar, Amu,
                                      nu, V, nu.e, ssq,
                                      R, keep, nprint, drawdelta,
                                      as.matrix(olddelta),  a, oldprob, ind, tau)
####################################################################

attributes(draws$taudraw)$class=c("bayesm.mat","mcmc")
attributes(draws$taudraw)$mcpar=c(1,R,keep)
if(drawdelta){
  attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
  attributes(draws$Deltadraw)$mcpar=c(1,R,keep)}
attributes(draws$betadraw)$class=c("bayesm.hcoef")
attributes(draws$nmix)$class="bayesm.nmix"

return(draws)
}