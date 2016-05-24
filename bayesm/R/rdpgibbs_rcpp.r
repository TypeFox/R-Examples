rDPGibbs=function(Prior,Data,Mcmc){
#
# Revision History: 
#   5/06 add rthetaDP
#   7/06 include rthetaDP in main body to avoid copy overhead
#   1/08 add scaling
#   2/08 add draw of lambda
#   3/08 changed nu prior support to dim(y) + exp(unif gird on nulim[1],nulim[2])
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
# purpose: do Gibbs sampling for density estimation using Dirichlet process model
#
# arguments:
#     Data is a list of y which is an n x k matrix of data
#     Prior is a list of (alpha,lambda,Prioralpha)
#       alpha: starting value
#       lambda_hyper: hyperparms of prior on lambda
#       Prioralpha: hyperparms of alpha prior; a list of (Istarmin,Istarmax,power)
#       if elements of the prior don't exist, defaults are assumed
#     Mcmc is a list of (R,keep,maxuniq)
#       R: number of draws
#       keep: thinning parameter
#       maxuniq: the maximum number of unique thetaStar values
#       nprint - print estimated time remaining on every nprint'th draw
#
# Output:
#     list with elements
#     alphadraw: vector of length R/keep, [i] is ith draw of alpha
#     Istardraw: vector of length R/keep, [i] is the number of unique theta's drawn from ith iteration
#     adraw
#     nudraw
#     vdraw
#     thetaNp1draws: list, [[i]] is ith draw of theta_{n+1}
#     inddraw: R x n matrix, [,i] is indicators of identity for each obs in ith iteration
#
# Model:
#        y_i ~ f(y|thetai)
#        thetai|G ~ G
#        G|lambda,alpha ~ DP(G|G0(lambda),alpha)
#
# Priors:
#        alpha: starting value
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
#
#

#  check arguments
#
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
alimdef=BayesmConstant.DPalimdef
nulimdef=BayesmConstant.DPnulimdef
vlimdef=BayesmConstant.DPvlimdef


if(missing(Prior)) {pandterm("requires Prior argument ")}
else
   {
    if(is.null(Prior$lambda_hyper)) {lambda_hyper=list(alim=alimdef,nulim=nulimdef,vlim=vlimdef)}
    else {lambda_hyper=Prior$lambda_hyper;
       if(is.null(lambda_hyper$alim)) {lambda_hyper$alim=alimdef}
       if(is.null(lambda_hyper$nulim)) {lambda_hyper$nulim=nulimdef} 
       if(is.null(lambda_hyper$vlim)) {lambda_hyper$vlim=vlimdef}
       }
    if(is.null(Prior$Prioralpha)) {Prioralpha=list(Istarmin=BayesmConstant.DPIstarmin,Istarmax=min(50,0.1*nobs),power=BayesmConstant.DPpower)}
    else {Prioralpha=Prior$Prioralpha;
       if(is.null(Prioralpha$Istarmin)) {Prioralpha$Istarmin=BayesmConstant.DPIstarmin} else {Prioralpha$Istarmin=Prioralpha$Istarmin}
       if(is.null(Prioralpha$Istarmax)) 
             {Prioralpha$Istarmax=min(50,0.1*nobs)} else {Prioralpha$Istarmax=Prioralpha$Istarmax}
       if(is.null(Prioralpha$power)) {Prioralpha$power=BayesmConstant.DPpower}
       }
   }
gamma= BayesmConstant.gamma
Prioralpha$alphamin=exp(digamma(Prioralpha$Istarmin)-log(gamma+log(nobs)))
Prioralpha$alphamax=exp(digamma(Prioralpha$Istarmax)-log(gamma+log(nobs)))
Prioralpha$n=nobs
#
# check Prior arguments for valdity
#
if(lambda_hyper$alim[1]<0) {pandterm("alim[1] must be >0")}
if(lambda_hyper$nulim[1]<0) {pandterm("nulim[1] must be >0")}
if(lambda_hyper$vlim[1]<0) {pandterm("vlim[1] must be >0")}
if(Prioralpha$Istarmin <1){pandterm("Prioralpha$Istarmin must be >= 1")}
if(Prioralpha$Istarmax <= Prioralpha$Istarmin){pandterm("Prioralpha$Istarmin must be > Prioralpha$Istarmax")}
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
    if(is.null(Mcmc$maxuniq)) {maxuniq=BayesmConstant.DPmaxuniq} else {maxuniq=Mcmc$maxuniq}
    if(is.null(Mcmc$SCALE)) {SCALE=BayesmConstant.DPSCALE} else {SCALE=Mcmc$SCALE}
    if(is.null(Mcmc$gridsize)) {gridsize=BayesmConstant.DPgridsize} else {gridsize=Mcmc$gridsize}
   }

#
# print out the problem
#
cat(" Starting Gibbs Sampler for Density Estimation Using Dirichlet Process Model",fill=TRUE)
cat(" ",nobs," observations on ",dimy," dimensional data",fill=TRUE)
cat(" ",fill=TRUE)
cat(" SCALE=",SCALE,fill=TRUE)
cat(" ",fill=TRUE)
cat(" Prior Parms: ",fill=TRUE)
cat("  G0 ~ N(mubar,Sigma (x) Amu^-1)",fill=TRUE)
cat("   mubar = ",0,fill=TRUE)
cat("   Sigma ~ IW(nu,nu*v*I)",fill=TRUE)
cat("   Amu ~ uniform[",lambda_hyper$alim[1],",",lambda_hyper$alim[2],"]",fill=TRUE)
cat("   nu ~ uniform on log grid on [",dimy-1+exp(lambda_hyper$nulim[1]),
             ",",dimy-1+exp(lambda_hyper$nulim[2]),"]",fill=TRUE)
cat("   v ~ uniform[",lambda_hyper$vlim[1],",",lambda_hyper$vlim[2],"]",fill=TRUE)
cat(" ",fill=TRUE)
cat("  alpha ~ (1-(alpha-alphamin)/(alphamax-alphamin))^power",fill=TRUE)
cat("   Istarmin = ",Prioralpha$Istarmin,fill=TRUE)
cat("   Istarmax = ",Prioralpha$Istarmax,fill=TRUE)
cat("   alphamin = ",Prioralpha$alphamin,fill=TRUE)
cat("   alphamax = ",Prioralpha$alphamax,fill=TRUE)
cat("   power = ",Prioralpha$power,fill=TRUE)
cat(" ",fill=TRUE)
cat(" Mcmc Parms: R= ",R," keep= ",keep," nprint= ",nprint," maxuniq= ",maxuniq," gridsize for lambda hyperparms= ",gridsize,
        fill=TRUE)
cat(" ",fill=TRUE)

###################################################################
# Wayne Taylor
# 1/29/2015
###################################################################
out = rDPGibbs_rcpp_loop(R,keep,nprint,
                         y, lambda_hyper, SCALE, maxuniq, Prioralpha, gridsize,
                         BayesmConstant.A,BayesmConstant.nuInc,BayesmConstant.DPalpha)
###################################################################

nmix=list(probdraw=matrix(c(rep(1,nrow(out$inddraw))),ncol=1),zdraw=out$inddraw,compdraw=out$thetaNp1draw)
attributes(nmix)$class="bayesm.nmix"
attributes(out$alphadraw)$class=c("bayesm.mat","mcmc")
attributes(out$Istardraw)$class=c("bayesm.mat","mcmc")
attributes(out$adraw)$class=c("bayesm.mat","mcmc")
attributes(out$nudraw)$class=c("bayesm.mat","mcmc")
attributes(out$vdraw)$class=c("bayesm.mat","mcmc")
return(list(alphadraw=out$alphadraw,Istardraw=out$Istardraw,adraw=out$adraw,nudraw=out$nudraw,
            vdraw=out$vdraw,nmix=nmix))
}
