rivDP = function(Data,Prior,Mcmc) {
#
# revision history:
#   P. Rossi 1/06
#   added draw of alpha 2/06
#   added automatic scaling 2/06
#   removed reqfun  7/07 -- now functions are in rthetaDP
#   fixed initialization of theta 3/09
#   fixed error in assigning user defined prior parms
#   W. Taylor 4/15 - added nprint option to MCMC argument
#
# purpose: 
#   draw from posterior for linear I.V. model with DP process for errors
#
# Arguments:
#   Data -- list of z,w,x,y
#        y is vector of obs on lhs var in structural equation
#        x is "endogenous" var in structural eqn
#        w is matrix of obs on "exogenous" vars in the structural eqn
#        z is matrix of obs on instruments
#   Prior -- list of md,Ad,mbg,Abg,mubar,Amu,nuV
#        md is prior mean of delta
#        Ad is prior prec
#        mbg is prior mean vector for beta,gamma
#        Abg is prior prec of same
#        lamda is a list of prior parms for DP draw
#              mubar is prior mean of means for "errors"
#              Amu is scale precision parm for means
#              nu,V parms for IW on Sigma (idential priors for each normal comp
#        alpha prior parm for DP process (weight on base measure)
#           or starting value if there is a prior on alpha (requires element Prioralpha)
#        Prioralpha list of hyperparms for draw of alpha (alphamin,alphamax,power,n)
#
#   Mcmc -- list of R,keep,starting values for delta,beta,gamma,theta
#        maxuniq is maximum number of unique theta values
#        R is number of draws
#        keep is thinning parameter
#        nprint - print estimated time remaining on every nprint'th draw
#        SCALE if scale data, def: TRUE
#        gridsize is the gridsize parm for alpha draws
#
#   Output: 
#      list of draws of delta,beta,gamma and thetaNp1 which is used for
#      predictive distribution of errors (density estimation)
# 
#   Model:
#
#    x=z'delta + e1
#    y=beta*x + w'gamma + e2
#        e1,e2 ~ N(theta_i)
#
#   Priors
#   delta ~ N(md,Ad^-1)
#   vec(beta,gamma) ~ N(mbg,Abg^-1)
#   theta ~ DPP(alpha|lambda)
#
#
#   extract data and check dimensios
#
if(missing(Data)) {pandterm("Requires Data argument -- list of z,w,x,y")}
    if(is.null(Data$w)) isgamma=FALSE else isgamma=TRUE
    if(isgamma) w = Data$w #matrix
    if(is.null(Data$z)) {pandterm("Requires Data element z")}
    z=Data$z
    if(is.null(Data$x)) {pandterm("Requires Data element x")}
    x=as.vector(Data$x)
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=as.vector(Data$y)
#
# check data for validity
#
n=length(y)
if(isgamma)
   {if(!is.matrix(w)) {pandterm("w is not a matrix")}
   dimg=ncol(w)
   if(n != nrow(w) ) {pandterm("length(y) ne nrow(w)")}}

if(!is.matrix(z)) {pandterm("z is not a matrix")}
dimd=ncol(z)
if(n != length(x) ) {pandterm("length(y) ne length(x)")}
if(n != nrow(z) ) {pandterm("length(y) ne nrow(z)")}


#
# extract elements corresponding to the prior
#
alimdef=BayesmConstant.DPalimdef
nulimdef=BayesmConstant.DPnulimdef
vlimdef=BayesmConstant.DPvlimdef

if(missing(Prior))
   {
    md=c(rep(0,dimd)) 
    Ad=diag(BayesmConstant.A,dimd) 
    if(isgamma) dimbg=1+dimg else dimbg=1
    mbg=c(rep(0,dimbg)) 
    Abg=diag(BayesmConstant.A,dimbg) 
 

    gamma= BayesmConstant.gamma  
    Istarmin=BayesmConstant.DPIstarmin
    alphamin=exp(digamma(Istarmin)-log(gamma+log(n)))
    Istarmax=floor(.1*n)
    alphamax=exp(digamma(Istarmax)-log(gamma+log(n)))
    power=BayesmConstant.DPpower
    Prioralpha=list(n=n,alphamin=alphamin,alphamax=alphamax,power=power)

    lambda_hyper=list(alim=alimdef,nulim=nulimdef,vlim=vlimdef)
   }

else  
   { 
    if(is.null(Prior$md)) md=c(rep(0,dimd)) else md=Prior$md
    if(is.null(Prior$Ad)) Ad=diag(BayesmConstant.A,dimd) else Ad=Prior$Ad
    if(isgamma) dimbg=1+dimg else dimbg=1
    if(is.null(Prior$mbg)) mbg=c(rep(0,dimbg)) else mbg=Prior$mbg
    if(is.null(Prior$Abg)) Abg=diag(BayesmConstant.A,dimbg) else Abg=Prior$Abg
    
    if(!is.null(Prior$Prioralpha))
       {Prioralpha=Prior$Prioralpha}
    else
       {gamma= BayesmConstant.gamma 
        Istarmin=BayesmConstant.DPIstarmin
        alphamin=exp(digamma(Istarmin)-log(gamma+log(n)))
        Istarmax=floor(.1*n)
        alphamax=exp(digamma(Istarmax)-log(gamma+log(n)))
        power=BayesmConstant.DPpower
        Prioralpha=list(n=n,alphamin=alphamin,alphamax=alphamax,power=power)}
    
    if(is.null(Prior$lambda_hyper)) {lambda_hyper=Prior$lambda_hyper}
    else
       {lambda_hyper=Prior$lambda_hyper;
          if(is.null(lambda_hyper$alim)) {lambda_hyper$alim=alimdef}
          if(is.null(lambda_hyper$nulim)) {lambda_hyper$nulim=nulimdef} 
          if(is.null(lambda_hyper$vlim)) {lambda_hyper$vlim=vlimdef}
       } 
   }
#
# check Prior arguments for valdity
#
if(lambda_hyper$alim[1]<0) {pandterm("alim[1] must be >0")}
if(lambda_hyper$nulim[1]<0) {pandterm("nulim[1] must be >0")}
if(lambda_hyper$vlim[1]<0) {pandterm("vlim[1] must be >0")}

#
# obtain starting values for MCMC
#
# we draw need inital values of delta, theta and indic
#

if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
theta=NULL
if(!is.null(Mcmc$delta)) 
   {delta = Mcmc$delta}
else
   {lmxz = lm(x~z,data.frame(x=x,z=z))
    delta = lmxz$coef[2:(ncol(z)+1)]}
if(!is.null(Mcmc$theta))
  {theta=Mcmc$theta }
else
  {onecomp=list(mu=c(0,0),rooti=diag(2))
   theta=vector("list",length(y))
   for(i in 1:n) {theta[[i]]=onecomp}
   }
dimd = length(delta)
if(is.null(Mcmc$maxuniq))
   {maxuniq=BayesmConstant.DPmaxuniq}
else
   {maxuniq=Mcmc$maxuniq}
if(is.null(Mcmc$R)) {pandterm("requres Mcmc argument, R")}
R = Mcmc$R
if(is.null(Mcmc$keep))
   {keep=BayesmConstant.keep}
else
   {keep=Mcmc$keep}
if(is.null(Mcmc$nprint))
{nprint=BayesmConstant.nprint}
else
{nprint=Mcmc$nprint}
if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
if(is.null(Mcmc$gridsize))
   {gridsize=BayesmConstant.DPgridsize}
else
   {gridsize=Mcmc$gridsize}
if(is.null(Mcmc$SCALE))
  {SCALE=BayesmConstant.DPSCALE}
else
  {SCALE=Mcmc$SCALE}

#
# scale and center
#
if(SCALE){
  scaley=sqrt(var(y))
  scalex=sqrt(var(x))
  meany=mean(y)
  meanx=mean(x)
  meanz=apply(z,2,mean)
  y=(y-meany)/scaley; x=(x-meanx)/scalex
  z=scale(z,center=TRUE,scale=FALSE)
  if(isgamma) {meanw=apply(w,2,mean);  w=scale(w,center=TRUE,scale=FALSE)}
}

#
# print out model
#
cat(" ",fill=TRUE)
cat("Starting Gibbs Sampler for Linear IV Model With DP Process Errors",fill=TRUE)
cat(" ",fill=TRUE)
cat(" nobs= ",n,"; ",ncol(z)," instruments",fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("mean of delta ",fill=TRUE)
print(md)
cat(" ",fill=TRUE)
cat("Adelta",fill=TRUE)
print(Ad)
cat(" ",fill=TRUE)
cat("mean of beta/gamma",fill=TRUE)
print(mbg)
cat(" ",fill=TRUE)
cat("Abeta/gamma",fill=TRUE)
print(Abg)
cat(" ",fill=TRUE)
cat("G0 ~ N(mubar,Sigma (x) Amu^-1)",fill=TRUE)
cat(" mubar = ",0,fill=TRUE)
cat(" Sigma ~ IW(nu,nu*v*I)",fill=TRUE)
cat(" Amu ~ uniform[",lambda_hyper$alim[1],",",lambda_hyper$alim[2],"]",fill=TRUE)
cat(" nu ~ uniform on log grid  [",2-1+exp(lambda_hyper$nulim[1]),
    ",",2-1+exp(lambda_hyper$nulim[2]),"]",fill=TRUE)
cat(" v ~ uniform[",lambda_hyper$vlim[1],",",lambda_hyper$vlim[2],"]",fill=TRUE)
cat("  ",fill=TRUE)
cat("Parameters of Prior on Dirichlet Process parm (alpha)",fill=TRUE)
cat("alphamin= ",Prioralpha$alphamin," alphamax= ",Prioralpha$alphamax," power=",
        Prioralpha$power,fill=TRUE)
cat("alpha values correspond to Istarmin = ",Istarmin," Istarmax = ",Istarmax,fill=TRUE)
cat(" ",fill=TRUE)
cat("MCMC parms: R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
cat("  maximum number of unique thetas= ",maxuniq,fill=TRUE)
cat("  gridsize for alpha draws= ",gridsize,fill=TRUE)
cat("  SCALE data= ",SCALE,fill=TRUE)
cat(" ",fill=TRUE)

###################################################################
# Wayne Taylor
# 3/14/2015
###################################################################
if(isgamma == FALSE) w=matrix()
out = rivDP_rcpp_loop(R,keep,nprint,dimd,mbg,Abg,md,Ad,y,isgamma,z,x,w, 
                      delta,PrioralphaList=Prioralpha,gridsize,SCALE,maxuniq,scalex,scaley,lambda_hyper,
                      BayesmConstant.A,BayesmConstant.nu)
###################################################################

nmix=list(probdraw=matrix(c(rep(1,length(out$thetaNp1draw))),ncol=1),zdraw=NULL,compdraw=out$thetaNp1draw)
#
# densitymix is in the format to be used with the generic mixture of normals plotting
# methods (plot.bayesm.nmix)
#
attributes(nmix)$class=c("bayesm.nmix")

attributes(out$deltadraw)$class=c("bayesm.mat","mcmc")
attributes(out$deltadraw)$mcpar=c(1,R,keep)
attributes(out$betadraw)$class=c("bayesm.mat","mcmc")
attributes(out$betadraw)$mcpar=c(1,R,keep)
attributes(out$alphadraw)$class=c("bayesm.mat","mcmc")
attributes(out$alphadraw)$mcpar=c(1,R,keep)
attributes(out$Istardraw)$class=c("bayesm.mat","mcmc")
attributes(out$Istardraw)$mcpar=c(1,R,keep)
if(isgamma){
  attributes(out$gammadraw)$class=c("bayesm.mat","mcmc")
  attributes(out$gammadraw)$mcpar=c(1,R,keep)}

if(isgamma) 
{ return(list(deltadraw=out$deltadraw,betadraw=out$betadraw,alphadraw=out$alphadraw,Istardraw=out$Istardraw,
              gammadraw=out$gammadraw,nmix=nmix))}
else
{ return(list(deltadraw=out$deltadraw,betadraw=out$betadraw,alphadraw=out$alphadraw,Istardraw=out$Istardraw,
              nmix=nmix))}
}