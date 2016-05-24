rnmixGibbs= function(Data,Prior,Mcmc){
#
# Revision History: 
#   P. Rossi 3/05
#   add check to see if Mubar is a vector  9/05
#   fixed bug in saving comps draw comps[[mkeep]]=  9/05
#   fixed so that ncomp can be =1; added check that nobs >= 2*ncomp   12/06
#   3/07 added classes
#   added log-likelihood  9/08
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
# purpose: do Gibbs sampling inference for a mixture of multivariate normals
#
# arguments:
#     Data is a list of y which is an n x k matrix of data -- each row
#       is an iid draw from the normal mixture
#     Prior is a list of (Mubar,A,nu,V,a,ncomp)
#       ncomp is required
#       if elements of the prior don't exist, defaults are assumed
#     Mcmc is a list of R, keep (thinning parameter), and nprint
# Output:
#     list with elements
#     pdraw -- R/keep x ncomp array of mixture prob draws
#     zdraw -- R/keep x nobs array of indicators of mixture comp identity for each obs
#     compsdraw -- list of R/keep lists of lists of comp parm draws
#        e.g. compsdraw[[i]] is ith draw -- list of ncomp lists
#             compsdraw[[i]][[j]] is list of parms for jth normal component
#             if jcomp=compsdraw[[i]][j]]
#                        ~N(jcomp[[1]],Sigma), Sigma = t(R)%*%R, R^{-1} = jcomp[[2]]
#
# Model:
#        y_i ~ N(mu_ind,Sigma_ind)
#        ind ~ iid multinomial(p)  p is a 1x ncomp vector of probs
# Priors:
#        mu_j ~ N(mubar,Sigma (x) A^-1)
#        mubar=vec(Mubar)
#        Sigma_j ~ IW(nu,V)
#        note: this is the natural conjugate prior -- a special case of multivariate 
#              regression
#        p ~ Dirchlet(a)
#
#  check arguments
#
#
# -----------------------------------------------------------------------------------------
llnmix=function(Y,z,comps){
  #
  # evaluate likelihood for mixture of normals
  #
  zu=unique(z)
  ll=0.0
  for(i in 1:length(zu)){
    Ysel=Y[z==zu[i],,drop=FALSE]
    ll=ll+sum(apply(Ysel,1,lndMvn,mu=comps[[zu[i]]]$mu,rooti=comps[[zu[i]]]$rooti))
  }
  return(ll)
}
# -----------------------------------------------------------------------------------------
if(missing(Data)) {pandterm("Requires Data argument -- list of y")}
if(is.null(Data$y)) {pandterm("Requires Data element y")}
y=Data$y
#
# check data for validity
#
if(!is.matrix(y)) {pandterm("y must be a matrix")}
nobs=nrow(y)
dimy=ncol(y)
#
# check for Prior
#
if(missing(Prior)) {pandterm("requires Prior argument ")}
else
{
  if(is.null(Prior$ncomp)) {pandterm("requires number of mix comps -- Prior$ncomp")}
  else {ncomp=Prior$ncomp}
  if(is.null(Prior$Mubar)) {Mubar=matrix(rep(0,dimy),nrow=1)} 
  else {Mubar=Prior$Mubar; if(is.vector(Mubar)) {Mubar=matrix(Mubar,nrow=1)}}
  if(is.null(Prior$A)) {A=matrix(BayesmConstant.A,ncol=1)} 
  else {A=Prior$A}
  if(is.null(Prior$nu)) {nu=dimy+BayesmConstant.nuInc} 
  else {nu=Prior$nu}
  if(is.null(Prior$V)) {V=nu*diag(dimy)} 
  else {V=Prior$V}
  if(is.null(Prior$a)) {a=c(rep(BayesmConstant.a,ncomp))}
  else {a=Prior$a}
}
#
# check for adequate no. of observations
#
if(nobs<2*ncomp)
{pandterm("too few obs, nobs should be >= 2*ncomp")}
#
# check dimensions of Priors
#
if(ncol(A) != nrow(A) || ncol(A) != 1)
{pandterm(paste("bad dimensions for A",dim(A)))}
if(!is.matrix(Mubar))
{pandterm("Mubar must be a matrix")}
if(nrow(Mubar) != 1 || ncol(Mubar) != dimy) 
{pandterm(paste("bad dimensions for Mubar",dim(Mubar)))}
if(ncol(V) != nrow(V) || ncol(V) != dimy)
{pandterm(paste("bad dimensions for V",dim(V)))}
if(length(a) != ncomp)
{pandterm(paste("a wrong length, length= ",length(a)))}
bada=FALSE
for(i in 1:ncomp){if(a[i] < 0) bada=TRUE}
if(bada) pandterm("invalid values in a vector")
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
  if(is.null(Mcmc$LogLike)) {LogLike=FALSE} else {LogLike=Mcmc$LogLike}
}

#
# print out the problem
#
cat(" Starting Gibbs Sampler for Mixture of Normals",fill=TRUE)
cat(" ",nobs," observations on ",dimy," dimensional data",fill=TRUE)
cat("     using ",ncomp," mixture components",fill=TRUE)
cat(" ",fill=TRUE)
cat(" Prior Parms: ",fill=TRUE)
cat("  mu_j ~ N(mubar,Sigma (x) A^-1)",fill=TRUE)
cat("  mubar = ",fill=TRUE)
print(Mubar)
cat("  precision parm for prior variance of mu vectors (A)= ",A,fill=TRUE)
cat("  Sigma_j ~ IW(nu,V) nu= ",nu,fill=TRUE)
cat("  V =",fill=TRUE)
print(V)
cat("  Dirichlet parameters ",fill=TRUE)
print(a)
cat(" ",fill=TRUE)
cat(" Mcmc Parms: R= ",R," keep= ",keep," nprint= ",nprint," LogLike= ",LogLike,fill=TRUE)

# pdraw=matrix(double(floor(R/keep)*ncomp),ncol=ncomp)
# zdraw=matrix(double(floor(R/keep)*nobs),ncol=nobs)
# compdraw=list()
compsd=list()
if(LogLike) ll=double(floor(R/keep))

#
# set initial values of z
#
z=rep(c(1:ncomp),(floor(nobs/ncomp)+1))
z=z[1:nobs]
cat(" ",fill=TRUE)
cat("starting value for z",fill=TRUE)
print(table(z))
cat(" ",fill=TRUE)
p=c(rep(1,ncomp))/ncomp # note this is not used
fsh()

#Wayne Taylor 8/18/14#####################################################
nmix = rnmixGibbs_rcpp_loop(y, Mubar, A, nu, V, a, p, z, R, keep, nprint);
##########################################################################

attributes(nmix)$class="bayesm.nmix"
  if(LogLike){
    zdraw = nmix$zdraw
    compdraw = nmix$compdraw
    ll = lapply(seq_along(compdraw), function(i) llnmix(y, zdraw[i,], compdraw[[i]]))
    return(list(ll=ll,nmix=nmix))
  }else{
    return(list(nmix=nmix))
  }
}