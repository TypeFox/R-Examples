# An example file fitting the FitzHugh-Nagumo equations to data in the
# new R Profiling Code. This will eventually be interfaced with the new R 
# "Partially Observed Markov Process (pomp)" class of objects. 

library(fda)
library(odesolve)
library(maxLik)
library(MASS)
library('Matrix')
library('SparseM')

source('../R/ProfileR.R')
source('../R/fhn.R')
source('exp2fhn.R')
source('../R/findif.ode.R')
source('../R/SSElik.R')
source('../R/SSEproc.R')
source('../R/makeid.R')
source('../R/genlin.R')
source('../R/cvar.R')
source('../R/Multinorm.R')
source('../R/Cproc.R')
source('../R/exp.Cproc.R')
source('../R/multinorm.shortcut.R')
source('../R/sse.shortcut.R')
source('../R/inneropt.R')
source('../R/makeexp.R')
source('../R/logtrans.R')

###############################
####   Data Generation  #######
###############################

t = seq(0,20,0.05)

pars = c(0.2,0.2,3)
names(pars) = c('a','b','c')
x0 = c(-1,1)
names(x0) = c('V','R')
fhn = make.fhn()
y = lsoda(x0,times=t,func=fhn$fn.ode,pars)
y = y[,2:3]

data = exp(y) + 0.05*array(rnorm(802),dim(y))


###############################
####  Basis Object      #######
###############################

knots = seq(0,20,0.2)
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	norder=norder,breaks=knots)


# Initial values for coefficients will be obtained by smoothing

fd.data = array(log(data),c(dim(data)[1],1,dim(data)[2]))

DEfd = data2fd(fd.data,t,bbasis,fdnames=list(NULL,NULL,c('V','R')) )

coefs = matrix(as.vector(DEfd$coefs),dim(DEfd$coefs)[1],dim(DEfd$coefs)[3])
colnames(coefs) = DEfd$fdnames[[3]]

###############################
####  Optimization Control  ###
###############################


control=list()                # Control parameters 
control$trace = 0
control$maxit = 1000
control$maxtry = 10
control$reltol = 1e-6
control$meth = "BFGS"

control.in = control
control.in$reltol = 1e-12
control.out = control
control.out$trace = 2

control.in$print.level = 0
control.in$iterlim = 1000



###############################
####     Optimization       ###
###############################

spars = c(0.2,0.2,2)          # Perturbed parameters
names(spars)=names(pars)
lambda = 10000

### SSE ####

res1 = Profile.sse(make.exp2fhn(),data,t,pars=spars,coefs,bbasis,lambda=lambda,in.meth='nlminb',out.meth='nls',
	control.in=control.in,control.out=control.out,pos=1,discrete=0)
	
### SSE with ProfileErr ###
profile.obj = sse.setup(pars=spars,coefs=coefs,fn=make.exp2fhn(),basisvals=bbasis,lambda=lambda,times=t,pos=1,discrete=0)
lik = profile.obj$lik
proc= profile.obj$proc

#res0 = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
#	control=control.out,times=t,data=data,lik=lik,proc=proc,pars=spars)
#
#ncoefs = array(res0$par,c(bbasis$nbasis,ncol(data)))
#
#res2 = optim(spars,ProfileErr,allpars=spars,times=t,data=data,coef=ncoefs,lik=lik,proc=proc,hessian=T,
#	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP,method="BFGS")
#
#res3 = nlminb(spars,ProfileErr,allpars=spars,times=t,data=data,coef=ncoefs,lik=lik,proc=proc,
#	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP)
#

Ires2 = inneropt(data,times=t,pars,coefs,lik,proc,in.meth='nlminb',control.in)

Ores2 = outeropt(data=data,times=t,pars=pars,coefs=coefs,lik=lik,proc=proc,in.meth="nlminb",out.meth="nlminb",
    control.in=control.in,control.out=control.out)




### Multinorm ####


var = c(1,0.01)
res4 = Profile.multinorm(make.exp2fhn(),data,t,pars=spars,coefs,bbasis,var=var, out.meth='BFGS', in.meth='nlminb',control.in=control.in,control.out=control.out,pos=1)

