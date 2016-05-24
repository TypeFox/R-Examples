#### Henon map examples
#
# This example demonstrates the use of CollocInfer on a discrete-time system.
# The Henon Map is a classical dynamical system that exhibits chaos.
# It's equations are given as
#
# x[t+1] = 1 - a*x[t]^2 + y[t]
# y[t+1] = b*x[t]

###############################
####   Data Generation  #######
###############################
hpars = c(1.4,0.3)

ntimes = 200

x = c(-1,1)
X = matrix(0,ntimes+20,2)
X[1,] = x

for(i in 2:(ntimes+20)){ X[i,] = make.Henon()$ode(i,X[i-1,],hpars,NULL) }

X = X[20+1:ntimes,]

plot(X,xlab='x',ylab='y',cex.lab=1.5,cex.axis=1.5)


Y = X + 0.05*matrix(rnorm(ntimes*2),ntimes,2)

t = 1:ntimes

coefs = as.matrix(Y)

###############################
####  Optimization Control  ###
###############################
control=list(trace = 0,maxit = 1000,maxtry = 10,reltol = 1e-6,meth = "BFGS")

control.in = control
control.in$reltol = 1e-12

control.out = control
control.out$trace = 2

###############################
####     Optimization       ###
###############################

hpars2 = c(1.3,0.4)          # Perturbed parameters
names(hpars2)=names(hpars)
lambda = 10000

### SSE for discrete process####

Ires1	= Smooth.LS(make.Henon(),data=Y,times=t,pars=hpars2,coefs=coefs,basisvals=NULL,
  lambda=lambda,in.meth='nlminb',control.in=control.in,discrete=TRUE)
  
Ores1 = Profile.LS(make.Henon(),data=Y,t,pars=hpars2,coefs,basisvals=NULL,
  lambda=lambda,in.meth='nlminb',out.meth='nls',control.in=control.in,
  control.out=control.out,discrete=TRUE)

### ProfileErr with SSEproc ####

profile.obj = LS.setup(pars=hpars2,coefs=coefs,fn=make.Henon(),basisvals=NULL,
  lambda=lambda,times=t,discrete=TRUE)
lik = profile.obj$lik
proc= profile.obj$proc
   
Ires2 = inneropt(data=Y,times=t,pars=hpars2,coefs,lik,proc,in.meth='nlminb',control.in)

Ores2 = outeropt(data=Y,times=t,pars=hpars2,coefs=coefs,lik=lik,proc=proc,
  in.meth="nlminb",out.meth="nlminb",control.in=control.in,control.out=control.out)

    
### Dproc ############

var = c(1,0.01)

Ires3 = Smooth.multinorm(make.Henon(),data=Y,t,pars=hpars2,coefs,basisvals=NULL,
  var=var,in.meth='nlminb',control.in=control.in,discrete=TRUE)

Ores3 = Profile.multinorm(fn=make.Henon(),data=Y,times=t,pars=hpars2,coefs=coefs,
  basisvals=NULL,var=var,fd.obj=NULL,more=NULL,quadrature=NULL,in.meth='nlminb',
  out.meth='optim',control.in=control.in,control.out=control.out,eps=1e-6,
  active=NULL,,discrete=TRUE)

