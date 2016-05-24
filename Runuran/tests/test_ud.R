#############################################################################
##                                                                         ##
##   Tests for special distributions                                       ##
##                                                                         ##
#############################################################################

## --- Load test routines and test parameters -------------------------------

source("test_routines.R")

## --- Auxiliary routines ---------------------------------------------------

## use u-error of method PINV for continuous distributions
cont.max.uerror.of <- function (distr, args) {
  uddist <- paste("ud",distr,"(",args,")", sep="")
  aqfunct <- eval(parse(text=paste("function(u) {uq(pinvd.new(",uddist,",uresolution=1e-12),u)}", sep="")))
  pfunct <- eval(parse(text=paste("function(x) {p",distr,"(x,",args,")}", sep="")))
  uerr <- unur.cont.uerror(samplesize,aqfunct,pfunct)
  if (uerr > 1.e-11) ## we do not want to be to stringent here
    stop(paste("max u-error exceeded for '",uddist,"': ",uerr,"\n",sep=""))
}

## use x-error of method DGT for discrete distributions
discr.max.xerror.of <- function (distr, args) {
  uddist <- paste("ud",distr,"(",args,",ub=10000)", sep="")
  aqfunct <- eval(parse(text=paste("function(u) {uq(dgtd.new(",uddist,"),u)}", sep="")))
  qfunct <- eval(parse(text=paste("function(x) {q",distr,"(x,",args,")}", sep="")))
  xerr <- unur.xerror(samplesize,aqfunct,qfunct)
  if (xerr > 1.e-10) ## we do not want to be to stringent here
    stop(paste("max x-error exceeded for '",uddist,"': ",xerr,"\n",sep=""))
}

## just run when neither CDF or quantile function is available 
cont.just.run <- function (distr, args) {
  uddist <- paste("ud",distr,"(",args,")", sep="")
  qfunct <- eval(parse(text=paste("function(u) { uq(pinvd.new(",uddist,",uresolution=1e-12),u) }", sep="")))
  x <- qfunct(runif(samplesize))
  if (! is.finite(sum(x)))
    stop(paste("inverse CDF for '",uddist,"': 'Inf' or 'NaN' occured\n",sep=""))
}

discr.just.run <- function (distr, args) {
  uddist <- paste("ud",distr,"(",args,",ub=10000)", sep="")
  qfunct <- eval(parse(text=paste("function(u) { uq(dgtd.new(",uddist,"),u) }", sep="")))
  x <- qfunct(runif(samplesize))
  if (! is.finite(sum(x)))
    stop(paste("inverse CDF for '",uddist,"': 'Inf' or 'NaN' occured\n",sep=""))
}


#############################################################################
##                                                                         ##
##  Test implementations of the ud<distr> functions.                       ##
##  The argments must coincide with the corresponding                      ##
##  d|p|q|r<distr> R functions.                                            ## 
##                                                                         ##
#############################################################################

## -- Continuous distributions ----------------------------------------------

cont.max.uerror.of("beta","shape1=3,shape2=7")
cont.max.uerror.of("cauchy","location=1,scale=2")
cont.just.run     ("chi","df=5")
cont.max.uerror.of("chisq","df=5")
cont.max.uerror.of("exp","rate=3")
cont.max.uerror.of("f","df1=5,df2=7")
cont.just.run     ("frechet","shape=3,location=2,scale=5")
cont.just.run     ("ghyp","lambda=1.5,alpha=3,beta=2,delta=1,mu=0")
cont.just.run     ("gig","theta=3,psi=1,chi=1")
cont.just.run     ("giga","theta=3,omega=2")
cont.just.run     ("giga","theta=3,omega=2,eta=2")
cont.just.run     ("gumbel","location=2,scale=5")
cont.just.run     ("hyperbolic","alpha=3,beta=2,delta=1,mu=0")
cont.just.run     ("ig","mu=3,lambda=2")
cont.just.run     ("laplace","location=2,scale=5")
cont.max.uerror.of("logis","location=2,scale=5")
cont.just.run     ("lomax","shape=2,scale=5")
cont.max.uerror.of("lnorm","mean=1,sd=2")
cont.just.run     ("meixner","alpha=0.03,beta=0.13,delta=0.57,mu=-0.001")
cont.max.uerror.of("norm","mean=1,sd=2")
cont.just.run     ("pareto","k=3,a=1")
## just.run     ("planck","a=3")  -- not implemented
cont.just.run     ("powerexp","shape=3")
cont.just.run     ("rayleigh","scale=3")
cont.just.run     ("slash","")
cont.max.uerror.of("t","df=2.5")
## max.uerror.of("triang","df=2.5") -- not implemented
cont.just.run     ("vg","lambda=2.25,alpha=210.,mu=0.001,beta=-5.1")
cont.just.run     ("vg","lambda=200.,alpha=210.,mu=0.001,beta=-5.1")
cont.max.uerror.of("weibull","shape=3,scale=2")

## -- Discrete distributions ----------------------------------------------

discr.max.xerror.of("binom","size=1000,prob=0.3456789")
discr.max.xerror.of("geom","prob=0.3456789")
discr.max.xerror.of("hyper","m=15,n=5,k=7")
discr.just.run("logarithmic","shape=0.3")
discr.max.xerror.of("nbinom","size=100,prob=0.3456789")
discr.max.xerror.of("pois","lambda=2.34567")

## --- End ------------------------------------------------------------------
