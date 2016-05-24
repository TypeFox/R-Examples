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
source('../R/findif.ode.R')
source('../R/SSElik.R')
source('../R/SSEproc.R')
source('../R/makeid.R')
source('../R/makeexp.R')
source('../R/genlin.R')
source('../R/cvar.R')
source('../R/Multinorm.R')
source('../R/Cproc.R')
source('../R/exp.Cproc.R')
source('../R/multinorm.shortcut.R')
source('../R/sse.shortcut.R')
source('../R/inneropt.R')

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

data = y + 0.05*array(rnorm(802),dim(y))



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

fd.data = array(data,c(dim(data)[1],1,dim(data)[2]))

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

#################################
### Initial Parameter Guesses ###
#################################

profile.obj = sse.setup(pars=pars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,lambda=10000,times=t)
lik = profile.obj$lik
proc= profile.obj$proc

pres = ParsMatchOpt(pars,coefs,proc)

npars = pres$pars

#############################################################
### If We Only Observe One State, We Can Re-Smooth Others ### 
#############################################################

tcoefs = coefs
tcoefs[,2] = 0

fres = FitMatchOpt(coefs=tcoefs,which=2,pars=pres$pars,proc)

ncoefs = fres$coefs

###############################
#### Parameter Optimization ###
###############################

spars = c(0.2,0.2,2)          # Perturbed parameters
names(spars)=names(pars)
lambda = 10000

### SSE Shortcuts ####

Ires1	= Smooth.sse(make.fhn(),data,t,pars=spars,coefs,bbasis,lambda=lambda,in.meth='nlminb',control.in=control.in)
  

Ores1 = Profile.sse(make.fhn(),data=data,times=t,pars=spars,coefs=coefs,basisvals=bbasis,lambda=lambda,in.meth='nlminb',out.meth='nls',
	control.in=control.in,control.out=control.out)
	

### SSE with ProfileErr ###


Ires2 = inneropt(data,times=t,pars,coefs,lik,proc,in.meth='nlminb',control.in)

Ores2 = outeropt(data=data,times=t,pars=pars,coefs=coefs,lik=lik,proc=proc,in.meth="nlminb",out.meth="nlminb",
    control.in=control.in,control.out=control.out)


### Multinorm ####

var = c(1,0.01)

Ires3 = Smooth.multinorm(make.fhn(),data,t,pars=spars,coefs,bbasis,var=var,in.meth='nlminb',control.in=control.in)

Ores3 = Profile.multinorm(make.fhn(),data,t,pars=spars,coefs,bbasis,var=var, 
  out.meth='nlminb', in.meth='nlminb',control.in=control.in,control.out=control.out)



# Lets look at the result

DEfd = fd(Ores1$coefs,bbasis)   # Data and reconstructed trajectory

par(mfrow=c(2,1))
 plotfit.fd(data,t,DEfd)
 

traj = as.matrix(Ores1$proc$bvals$bvals%*%Ores1$coefs)    # Look at how well the
colnames(traj) = Ores1$proc$more$names                    # derivative of the
dtraj = as.matrix(Ores1$proc$bvals$dbvals%*%Ores1$coefs)  # trajectory fits the 
ftraj = Ores1$proc$more$fn(t,traj,Ores1$pars)             # right hand side. 

matplot(dtraj,type='l',col=2)
matplot(ftraj,type='l',col=4,add=TRUE) 
 
 
Profile.covariance(pars=Ores1$pars,times=t,data=data,coefs=Ores1$coefs,
  lik=Ores1$lik,proc=Ores1$proc)
  
  
  