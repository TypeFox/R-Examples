# An example file fitting the FitzHugh-Nagumo equations to data in the
# new R Profiling Code. This will eventually be interfaced with the new R 
# "Partially Observed Markov Process (pomp)" class of objects. 

# The equations are a simplification of a famous set of equations developed
# by Hodgkin and Huxley in 1952 to describe the time course of the action
# potential in a squid neuron.

# The variables are:
#  V ... voltage or potential difference across the membrane of the neurone
#  R ... a combination of the effects of three ion channels that restore
#        the resting potential holding between spikes.
#  V is potentially measureable, but R is not since it simplifes ther
#    behavior of the actual ion channels.

#  The equations are;
#  DV = c(V - V^3/3 + R)
#  DR = -(b R + V - a)/c

#  The R equation is a simple first order linear constant coefficient dynamic
#  system with rate b/c and forced by (V - a)/c
#  The V equation is nonlinear due to the -c V^3 term, and is forced by c R.

#  The data are simulated by adding Gaussian error to a solution of the
#  differential equation.

#  Last modifed 2 Feb 2010 by Jim

library("CollocInfer")

#  Set up 41 time values rangeing from 0 to 20 in steps of 0.5

t = seq(0,20,0.5)

n = length(t)

#  Set up parameter values; a = b = 0.2 and c = 3.0

pars = c(0.2,0.2,3)
names(pars) = c('a','b','c')

#  initial values and variable labels

x0 = c(-1,1)
names(x0) = c('V','R')

#  Define the functions that evaluate the right sides of the equations
#  as well as various derivatives

fhn = make.fhn()

#  approximate the solution for these values over about three cycles
#  Here we use the function lsoda found in package odesolve that 
#  approximates the solution of the differential equation given a pair
#  of initial values.

solution = lsoda(x0,times=t,func=fhn$fn.ode,pars)

#  extract the function values at time points

y = solution[,2:3]

#  add noise to the function values to define data

sigerr = 0.5
data = y + sigerr*array(rnorm(2*n),dim(y))

#  put data for R variable to NA to indicate that R is not measured

data[,2] = NA

#  plot the solution for both variables and add points for V

matplot(t, y, "l")
# matpoints(t, data, pch="*")
points(t, data[,1], pch="*")

###############################
####  Basis object bbasis  ####
###############################

knots  = seq(0,20,0.2)
norder = 4
nbasis = length(knots) + norder - 2
range  = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	                            norder=norder,breaks=knots)

########################################
### Initial values for coefficients  ###
########################################

fd.data = array(data,c(n,1,2))

bfdPar = fdPar(bbasis, 2, 1e-2)

Vfd = smooth.basis(t,fd.data[,,1],bfdPar)$fd
plotfit.fd(fd.data[,,1], t, Vfd)

coefs0 = matrix(0,nbasis,2)
coefs0[,1] = Vfd$coefs

DEfd = fd(coefs0, bbasis, fdnames=list(NULL,NULL,c('V','R')))

######################################################################
### Set up likelihood and process functions for squared error loss ###
######################################################################

#  set up functions to evaluate values of ODE right side and values
#  of its derivatives

profile.obj = LS.setup(pars=pars, fn=fhn, lambda=1000, times=t,
                       fd.obj=DEfd)

lik  = profile.obj$lik

proc = profile.obj$proc

#############################################
####  set optimization control parameters ###
#############################################

#  control

control=list()                 

#  control.in

control.in = control
control.in$trace    = 0
control.in$rel.tol  = 1e-4
control.in$iter.max = 100

#  control.out

control.out = control
control.out$trace = 1
control.out$rel.tol = 1e-4

#####################################################
#  Estimate coefficients for R by fitting equation  #
#####################################################

res1 = FitMatchOpt(coefs0, which=2, pars, proc ,meth='nlminb',
                   control=control.in)
coefs1 = res1$coefs

#  define current fit to equation by the coefficients returned by 
#  FitMatchOpt

DEfd$coefs = coefs1

plot(DEfd)
lines(t, y[,2], lty=2)

###############################
#### Parameter Optimization ###
###############################

# Perturbed parameter starting values

spars = pars*exp(rnorm(3)*0.2)          
names(spars) = names(pars)

#  set smoothing level (value to be used for both variables)

lambda = 1

####################################################
### Estimate parameters using squared error loss ###
####################################################

control.in$trace    = 1 
Ires1 = inneropt(data=data,times=t,pars=spars,coefs=coefs1,lik,proc,
                 in.meth='nlminb',control.in=control.in)

control.in$trace    = 0
Ores1 = Profile.LS(fhn,data=data,times=t,pars=spars,coefs=coefs1,
                    basisvals=bbasis,lambda=lambda,names=names(x0),
                    in.meth='nlminb',out.meth='optim',
	              control.in=control.in,control.out=control.out)
	
####################################################
###         SSE with ProfileErr                  ###
####################################################

control.in$trace    = 0
control.out$trace   = 6

Ores2 = outeropt(data=data,times=t,pars=spars,coefs=coefs1,
                 lik=lik,proc=proc,
                 in.meth="nlminb",out.meth="optim",
                 control.in=control.in,control.out=control.out)

####################################################
###          Lets look at the result             ###
####################################################

DEfd = fd(Ores2$coefs,bbasis)   # Data and reconstructed trajectory

#  plot estimated and true solutions to the equation

plot(DEfd)
matlines(t,y,lty=2)

#  plot fit to V data

plotfit.fd(data[,1],t,DEfd[1])
lines(t,y[,1],lty=2)
 
#  Ores2 does not return lik and proc, 
#  so Ores$lik and Ores$proc changed to lik and proc, resp.

traj = as.matrix(proc$bvals$bvals   %*% Ores2$coefs) # Look at how well the
colnames(traj) = proc$more$names                     # derivative of the
dtraj = as.matrix(proc$bvals$dbvals %*% Ores2$coefs) # trajectory fits the 
ftraj = proc$more$fn(t,traj,Ores2$pars)              # right hand side. 

matplot(dtraj,type='l',col=2)
matplot(ftraj,type='l',col=4,add=TRUE) 
 
#  This doesn't work

Profile.covariance(pars=Ores$pars,times=t,data=data,coefs=Ores$coefs,lik=lik,proc=proc)
  
  
	
