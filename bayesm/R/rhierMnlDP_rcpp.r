rhierMnlDP=function(Data,Prior,Mcmc){
#
#  created 3/08 by Rossi from rhierMnlRwMixture adding DP draw for to replace finite mixture of normals
#
# revision history:
#   changed 12/17/04 by rossi to fix bug in drawdelta when there is zero/one unit
#   in a mixture component
#   added loglike output, changed to reflect new argument order in llmnl, mnlHess 9/05
#   changed weighting scheme to (1-w)logl_i + w*Lbar (normalized) 12/05
#   3/07 added classes
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
# purpose: run hierarchical mnl logit model with mixture of normals 
#   using RW and cov(RW inc) = (hess_i + Vbeta^-1)^-1
#   uses normal approximation to pooled likelihood
#
# Arguments:
#   Data contains a list of (p,lgtdata, and possibly Z)
#      p is number of choice alternatives
#      lgtdata is a list of lists (one list per unit)
#          lgtdata[[i]]=list(y,X)
#             y is a vector indicating alternative chosen
#               integers 1:p indicate alternative
#             X is a length(y)*p x nvar matrix of values of
#               X vars including intercepts
#             Z is an length(lgtdata) x nz matrix of values of variables
#               note: Z should NOT contain an intercept
#   Prior contains a list of (deltabar,Ad,lambda_hyper,Prioralpha)
#       alpha: starting value
#       lambda_hyper: hyperparms of prior on lambda
#       Prioralpha: hyperparms of alpha prior; a list of (Istarmin,Istarmax,power)
#       if elements of the prior don't exist, defaults are assumed
#   Mcmc contains a list of (s,c,R,keep,nprint)
#
# Output:  as list containing
#   Deltadraw R/keep  x nz*nvar matrix of draws of Delta, first row is initial value
#   betadraw is nlgt x nvar x R/keep array of draws of betas
#   probdraw is R/keep x 1 matrix of draws of probs of mixture components
#   compdraw is a list of list of lists (length R/keep)
#      compdraw[[rep]] is the repth draw of components for mixtures
#   loglike  log-likelikelhood at each kept draw
#
# Priors:
#    beta_i = D %*% z[i,] + u_i
#       vec(D)~N(deltabar)
#       u_i ~ N(theta_i)
#       theta_i~G
#       G|lambda,alpha ~ DP(G|G0(lambda),alpha)
#
#        lambda:
#           G0 ~ N(mubar,Sigma (x) Amu^-1)
#           mubar=vec(mubar)
#           Sigma ~ IW(nu,nu*v*I)  note: mode(Sigma)=nu/(nu+2)*v*I
#           mubar=0
#           amu is uniform on grid specified by alim
#           nu is log uniform, nu=d-1+exp(Z) z is uniform on seq defined bvy nulim
#           v is uniform on sequence specificd by vlim
#
#        Prioralpha:
#           alpha ~ (1-(alpha-alphamin)/(alphamax-alphamin))^power
#           alphamin=exp(digamma(Istarmin)-log(gamma+log(N)))
#           alphamax=exp(digamma(Istarmax)-log(gamma+log(N)))
#           gamma= .5772156649015328606
#
# MCMC parameters
#   s is the scaling parameter for the RW inc covariance matrix; s^2 Var is inc cov
#      matrix
#   w is parameter for weighting function in fractional likelihood
#      w is the weight on the normalized pooled likelihood 
#   R is number of draws
#   keep is thinning parameter, keep every keepth draw
#   nprint - print estimated time remaining on every nprint'th draw
#--------------------------------------------------------------------------------------------------

llmnlFract=
function(beta,y,X,betapooled,rootH,w,wgt){
z=as.vector(rootH%*%(beta-betapooled))
return((1-w)*llmnl(beta,y,X)+w*wgt*(-.5*(z%*%z)))
}

#
#  check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")}
  if(is.null(Data$p)) {pandterm("Requires Data element p (# chce alternatives)") }
  p=Data$p
  if(is.null(Data$lgtdata)) {pandterm("Requires Data element lgtdata (list of data for each unit)")}
  lgtdata=Data$lgtdata
  nlgt=length(lgtdata)
  drawdelta=TRUE
if(is.null(Data$Z)) { cat("Z not specified",fill=TRUE); fsh() ; drawdelta=FALSE}
  else {if (nrow(Data$Z) != nlgt) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number logits ",nlgt))}
      else {Z=Data$Z}}
  if(drawdelta) {
     nz=ncol(Z)
     colmeans=apply(Z,2,mean)
     if(sum(colmeans) > .00001) 
       {pandterm(paste("Z does not appear to be de-meaned: colmeans= ",colmeans))}
  }
#
# check lgtdata for validity
#
ypooled=NULL
Xpooled=NULL
if(!is.null(lgtdata[[1]]$X)) {oldncol=ncol(lgtdata[[1]]$X)}
for (i in 1:nlgt) 
{
    if(is.null(lgtdata[[i]]$y)) {pandterm(paste("Requires element y of lgtdata[[",i,"]]"))}
    if(is.null(lgtdata[[i]]$X)) {pandterm(paste("Requires element X of lgtdata[[",i,"]]"))}
    ypooled=c(ypooled,lgtdata[[i]]$y)
    nrowX=nrow(lgtdata[[i]]$X)
    if((nrowX/p) !=length(lgtdata[[i]]$y)) {pandterm(paste("nrow(X) ne p*length(yi); exception at unit",i))}
    newncol=ncol(lgtdata[[i]]$X)
    if(newncol != oldncol) {pandterm(paste("All X elements must have same # of cols; exception at unit",i))}
    Xpooled=rbind(Xpooled,lgtdata[[i]]$X)
    oldncol=newncol
}
nvar=ncol(Xpooled)
levely=as.numeric(levels(as.factor(ypooled)))
if(length(levely) != p) {pandterm(paste("y takes on ",length(levely)," values -- must be = p"))}
bady=FALSE
for (i in 1:p )
{
    if(levely[i] != i) bady=TRUE
}
cat("Table of Y values pooled over all units",fill=TRUE)
print(table(ypooled))
if (bady) 
  {pandterm("Invalid Y")}
#
# check on prior
#
alimdef=BayesmConstant.DPalimdef
nulimdef=BayesmConstant.DPnulimdef
vlimdef=BayesmConstant.DPvlimdef

if(missing(Prior)) {Prior=NULL}

if(is.null(Prior$lambda_hyper)) {lambda_hyper=list(alim=alimdef,nulim=nulimdef,vlim=vlimdef)}
   else {lambda_hyper=Prior$lambda_hyper;
       if(is.null(lambda_hyper$alim)) {lambda_hyper$alim=alimdef}
       if(is.null(lambda_hyper$nulim)) {lambda_hyper$nulim=nulimdef} 
       if(is.null(lambda_hyper$vlim)) {lambda_hyper$vlim=vlimdef}
       }
if(is.null(Prior$Prioralpha)) {Prioralpha=list(Istarmin=BayesmConstant.DPIstarmin,Istarmax=min(50,0.1*nlgt),power=BayesmConstant.DPpower)}
   else {Prioralpha=Prior$Prioralpha;
       if(is.null(Prioralpha$Istarmin)) {Prioralpha$Istarmin=BayesmConstant.DPIstarmin} else {Prioralpha$Istarmin=Prioralpha$Istarmin}
       if(is.null(Prioralpha$Istarmax)) 
       {Prioralpha$Istarmax=min(50,0.1*nlgt)} else {Prioralpha$Istarmax=Prioralpha$Istarmax}
      if(is.null(Prioralpha$power)) {Prioralpha$power=BayesmConstant.DPpower}
   }	
gamma= BayesmConstant.gamma
Prioralpha$alphamin=exp(digamma(Prioralpha$Istarmin)-log(gamma+log(nlgt)))
Prioralpha$alphamax=exp(digamma(Prioralpha$Istarmax)-log(gamma+log(nlgt)))
Prioralpha$n=nlgt
#
# check Prior arguments for valdity
#
if(lambda_hyper$alim[1]<0) {pandterm("alim[1] must be >0")}
if(lambda_hyper$nulim[1]<0) {pandterm("nulim[1] must be >0")}
if(lambda_hyper$vlim[1]<0) {pandterm("vlim[1] must be >0")}
if(Prioralpha$Istarmin <1){pandterm("Prioralpha$Istarmin must be >= 1")}
if(Prioralpha$Istarmax <= Prioralpha$Istarmin){pandterm("Prioralpha$Istarmin must be < Prioralpha$Istarmax")}	

if(is.null(Prior$Ad) & drawdelta) {Ad=BayesmConstant.A*diag(nvar*nz)} else {Ad=Prior$Ad}
if(drawdelta) {if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}}
if(is.null(Prior$deltabar)& drawdelta) {deltabar=rep(0,nz*nvar)} else {deltabar=Prior$deltabar}
  if(drawdelta) {if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}}
#
# check on Mcmc
#
if(missing(Mcmc)) 
  {pandterm("Requires Mcmc list argument")}
else 
   { 
    if(is.null(Mcmc$s)) {s=BayesmConstant.RRScaling/sqrt(nvar)} else {s=Mcmc$s}
    if(is.null(Mcmc$w)) {w=BayesmConstant.w}  else {w=Mcmc$w}
    if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
    if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
      if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')} 
    if(is.null(Mcmc$maxuniq)) {maxuniq=BayesmConstant.DPmaxuniq} else {keep=Mcmc$maxuniq}
    if(is.null(Mcmc$gridsize)) {gridsize=BayesmConstant.DPgridsize} else {gridsize=Mcmc$gridsize}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {R=Mcmc$R}
    }
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting MCMC Inference for Hierarchical Logit:",fill=TRUE)
cat("   Dirichlet Process Prior",fill=TRUE)
cat(paste("  ",p," alternatives; ",nvar," variables in X"),fill=TRUE)
cat(paste("   for ",nlgt," cross-sectional units"),fill=TRUE)
cat(" ",fill=TRUE)
cat(" Prior Parms: ",fill=TRUE)
cat("  G0 ~ N(mubar,Sigma (x) Amu^-1)",fill=TRUE)
cat("   mubar = ",0,fill=TRUE)
cat("   Sigma ~ IW(nu,nu*v*I)",fill=TRUE)
cat("   Amu ~ uniform[",lambda_hyper$alim[1],",",lambda_hyper$alim[2],"]",fill=TRUE)
cat("   nu ~ uniform on log grid  [",nvar-1+exp(lambda_hyper$nulim[1]),
             ",",nvar-1+exp(lambda_hyper$nulim[2]),"]",fill=TRUE)
cat("   v ~ uniform[",lambda_hyper$vlim[1],",",lambda_hyper$vlim[2],"]",fill=TRUE)
cat(" ",fill=TRUE)
cat("  alpha ~ (1-(alpha-alphamin)/(alphamax-alphamin))^power",fill=TRUE)
cat("   Istarmin = ",Prioralpha$Istarmin,fill=TRUE)
cat("   Istarmax = ",Prioralpha$Istarmax,fill=TRUE)
cat("   alphamin = ",Prioralpha$alphamin,fill=TRUE)
cat("   alphamax = ",Prioralpha$alphamax,fill=TRUE)
cat("   power = ",Prioralpha$power,fill=TRUE)
cat(" ",fill=TRUE)
if(drawdelta) 
{
   cat("deltabar",fill=TRUE)
   print(deltabar)
   cat("Ad",fill=TRUE)
   print(Ad)
}
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat(paste("s=",round(s,3)," w= ",w," R= ",R," keep= ",keep," nprint= ",nprint," maxuniq= ",maxuniq,
          " gridsize for lambda hyperparms= ",gridsize),fill=TRUE)
cat("",fill=TRUE)
#
# allocate space for draws
#
oldbetas=matrix(double(nlgt*nvar),ncol=nvar)

#
# intialize compute quantities for Metropolis
#
cat("initializing Metropolis candidate densities for ",nlgt," units ...",fill=TRUE)
fsh()
#
#  now go thru and computed fraction likelihood estimates and hessians
#
#       Lbar=log(pooled likelihood^(n_i/N))
#
#       fraction loglike = (1-w)*loglike_i + w*Lbar
#
betainit=c(rep(0,nvar))
#
#  compute pooled optimum
#
out=optim(betainit,llmnl,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-6), 
     X=Xpooled,y=ypooled)
betapooled=out$par
H=mnlHess(betapooled,ypooled,Xpooled)
rootH=chol(H)
#
# initialize betas for all units
#
for (i in 1:nlgt) 
{
   wgt=length(lgtdata[[i]]$y)/length(ypooled)
   out=optim(betapooled,llmnlFract,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-4), 
   X=lgtdata[[i]]$X,y=lgtdata[[i]]$y,betapooled=betapooled,rootH=rootH,w=w,wgt=wgt)
   if(out$convergence == 0) 
     { hess=mnlHess(out$par,lgtdata[[i]]$y,lgtdata[[i]]$X)
       lgtdata[[i]]=c(lgtdata[[i]],list(converge=1,betafmle=out$par,hess=hess)) }
   else
     { lgtdata[[i]]=c(lgtdata[[i]],list(converge=0,betafmle=c(rep(0,nvar)),
        hess=diag(nvar))) }
   oldbetas[i,]=lgtdata[[i]]$betafmle
   if(i%%50 ==0) cat("  completed unit #",i,fill=TRUE)
   fsh()
}

#Initialize placeholders when drawdelta == FALSE
if (drawdelta==FALSE){
  Z = matrix(0)
  deltabar = 0
  Ad = matrix(0)
}

###################################################################
# Wayne Taylor
# 2/21/2015
###################################################################
out = rhierMnlDP_rcpp_loop(R,keep,nprint,
                           lgtdata,Z,deltabar,Ad,Prioralpha,lambda_hyper,
                           drawdelta,nvar,oldbetas,s,maxuniq,gridsize,
                           BayesmConstant.A,BayesmConstant.nuInc,BayesmConstant.DPalpha)
###################################################################

if(drawdelta){
  attributes(out$Deltadraw)$class=c("bayesm.mat","mcmc")
  attributes(out$Deltadraw)$mcpar=c(1,R,keep)}
attributes(out$betadraw)$class=c("bayesm.hcoef")
attributes(out$nmix)$class="bayesm.nmix"
attributes(out$adraw)$class=c("bayesm.mat","mcmc")
attributes(out$nudraw)$class=c("bayesm.mat","mcmc")
attributes(out$vdraw)$class=c("bayesm.mat","mcmc")
attributes(out$Istardraw)$class=c("bayesm.mat","mcmc")
attributes(out$alphadraw)$class=c("bayesm.mat","mcmc")

return(out)
}