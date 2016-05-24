##########################################
##       FitzHugh-Nagumo Example       ###
##########################################

# Obtain some pre-generated data

FhNdata

# Now we want to set up a basis expansion; this makes use of
# routines in the fda library

knots  = seq(0,20,0.5)
norder = 3
nbasis = length(knots) + norder - 2
range  = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	norder=norder,breaks=knots)
	
	
# We'll start off by creating a smooth of each state variable to get initial
# values for the parameters.

bfdPar = fdPar(bbasis,lambda=1,int2Lfd(1))
DEfd = smooth.basis(FhNtimes,FhNdata,bfdPar,fdnames=list(NULL,c('V','R')) )$fd

coefs = DEfd$coefs
colnames(coefs) = FhNvarnames

#### Option 1: pre-defined set-up for squared error

lambda = 1000

# First just optimize the coefficients

Ires1	= Smooth.LS(make.fhn(),FhNdata,FhNtimes,FhNpars,coefs,bbasis,lambda,
  in.meth='nlminb')

# Let's have a look at this

IDEfd1 = fd(Ires1$coefs,bbasis)
plot(IDEfd1)
matplot(FhNtimes,FhNdata,add=TRUE)


# Now we can do the profiling

Ores1 = Profile.LS(make.fhn(),FhNdata,FhNtimes,FhNpars,coefs=coefs,
  basisvals=bbasis,lambda=lambda)

# And look at the result

ODEfd1 = fd(Ores1$coefs,bbasis)
plot(ODEfd1)
matplot(FhNtimes,FhNdata,add=TRUE)
  
# If we only coded the FitzHugh Nagumo funciton with no derivatives we would just
# have make.fhn()$fn, in which case we default to estimating the derivatives by
# finite differencing

Ires1a	= Smooth.LS(make.fhn()$fn,FhNdata,FhNtimes,FhNpars,coefs,bbasis,lambda,
  in.meth='nlminb')

# Now we can do the profiling

Ores1a = Profile.LS(make.fhn()$fn,FhNdata,FhNtimes,FhNpars,coefs=coefs,
  basisvals=bbasis,lambda=lambda)
  
  
  
#### Option 2:  set-up functions for proc and lik objects and then call
# optimizers

profile.obj = LS.setup(pars=FhNpars,fn=make.fhn(),lambda=lambda,times=FhNtimes,
      coefs=coefs,basisvals=bbasis)
lik = profile.obj$lik
proc= profile.obj$proc

# Now we can get initial parameter estimates from "gradient matching"

pres = ParsMatchOpt(FhNpars,coefs,proc)
npars = pres$pars

# Smoothing can be done more specifically with 'inneropt'

Ires2 = inneropt(FhNdata,times=FhNtimes,npars,coefs,lik,proc)

# And we can also do profiling

Ores2 = outeropt(FhNdata,FhNtimes,npars,coefs,lik,proc)

#### Option 3:  set everything up manually

## lik object

lik2 = make.SSElik()                      # we are using squared error
lik2$bvals = eval.basis(FhNtimes,bbasis)  # values of the basis at observation times
lik2$more = make.id()                     # identity transformation
lik2$more$weights = c(1,0)                # only use data for V

## proc object

qpts = knots[1:(length(knots)-1)]+diff(knots)/2  # Collocation points at midpoints
                                                 # between knots

proc2 = make.SSEproc()                    # we are using squared error
proc2$bvals = list(bvals = eval.basis(qpts,bbasis),    # values and derivative of
                  dbvals = eval.basis(qpts,bbasis,1))  # basis at collocation points
proc2$more = make.fhn()                   # FitzHugh-Nagumo right hand side
proc2$more$names = FhNvarnames            # State variable names
proc2$more$parnames = FhNparnames         # Parameter names
proc2$more$qpts = qpts                    # Collocation points
proc2$more$weights  = c(1000,1000)        # Weight relative to observations of
                                          # each collocation point.
                                          
# We can also specify optimization methods

in.meth = 'nlminb'            # Inner Optimization
control.in = list(rel.tol=1e-12,iter.max=1000,eval.max=2000,trace=0)  # Optimization control

out.meth = 'optim'            # Outer Optimization
control.out = list(trace=6,maxit=100,reltol=1e-8,meth='BFGS') # Optimization control
                                          
                                          
# Now we will also remove the 'R' part of the data (ie, assume we did not measure
# it) and set the corresponding coefficients to zero

new.FhNdata = FhNdata
new.FhNdata[,2] = NA

new.coefs = coefs
new.coefs[,2] = 0

# We can now fix the smooth or 'V' and try and choose 'R' to make the differential
# equation match as well as possible.

fres = FitMatchOpt(coefs=new.coefs,which=2,pars=FhNpars,proc2)
new.coefs2 = fres$coefs

# And we can call the same inner optimization as above, this time with our
# own control parameters and method

Ires3 = inneropt(new.FhNdata,FhNtimes,FhNpars,new.coefs2,lik2,proc2,
  in.meth=in.meth,control.in=control.in)

# And we can also do profiling with specified controls

Ores3 = outeropt(new.FhNdata,FhNtimes,FhNpars,coefs=new.coefs2,lik=lik2,proc=proc2,
  in.meth=in.meth,out.meth=out.meth,control.in=control.in,control.out=control.out)

# We can also use finite differencing to calculate derivatives of the right hand side of the
# FitzHugh-Nagumo equations, this is defined by modifying the proc object.

proc2a = make.SSEproc()                    # we are using squared error
proc2a$bvals = list(bvals = eval.basis(qpts,bbasis),    # values and derivative of
                  dbvals = eval.basis(qpts,bbasis,1))  # basis at collocation points
proc2a$more = make.findif.ode()                   # FitzHugh-Nagumo right hand side
proc2a$more$names = FhNvarnames            # State variable names
proc2a$more$parnames = FhNparnames         # Parameter names
proc2a$more$qpts = qpts                    # Collocation points
proc2a$more$weights  = c(1000,1000)        # Weight relative to observations of
                                           # each collocation point.
proc2a$more$more = list(fn = make.fhn()$fn, eps=1e-6) # Tell findif that the function
                                                      # to difference is fhn and the
                                                      # difference stepsize


Ires3a = inneropt(new.FhNdata,FhNtimes,FhNpars,new.coefs2,lik2,proc2a,
  in.meth=in.meth,control.in=control.in)

# And we can also do profiling with specified controls

Ores3a = outeropt(new.FhNdata,FhNtimes,FhNpars,coefs=new.coefs2,lik=lik2,proc=proc2a,
  in.meth=in.meth,out.meth=out.meth,control.in=control.in,control.out=control.out)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##########
# This demonstration illustrates the use of log transformations with
# CollocInfer. We use the SEIR equations with a seasonally varying
# infection rate for this purpose.
#
# The equations are given as
#
# Sdot = mu - beta(t)*S*(I+i) - nu*S        (Susceptibles)
# Edot = beta(t)*S*(I+i) - (sigma+nu)*E     (Exposed)
# Idot = sigma*E - (gamma+nu)*I              (Infectious)
# 
# Here beta(t) - the infection rate - is parameterized by a sinusoidal function
# plus a constant.
#
# Traditionally, there is an additional state
#
# Rdot = gamma*I - nu*R
#
# However, since we only observe I, R contributes nothing to the data fit and
# we have removed it from the system. 
#
# Other parameters are
# i - a visiting process
# nu - death rate
# sigma - the rate of movement from Exposed to Infectious.
# gamma - the rate of recovery from infection.
#
# It is generally more stable to solve the equations for the log states rather
# than the states themselves. CollocInfer contains a number of useful tools
# that let you make this transition without needing to re-code all of your
# differential equations.

library('CollocInfer')
            
            
#### Get some data and parameters

SEIRvarnames = SEIRvarnames
SEIRparnames = SEIRparnames

SEIRtimes = SEIRtimes
SEIRdata = SEIRdata

SEIRpars = SEIRpars


#### Now format the data so that S and E measurements are listed as NA

data = cbind(matrix(NA,length(SEIRdata),2),SEIRdata)

# We'll also look at the log observations

logdata = log(data)



#### define the right side evaluation function

SEIRfn = make.SEIR()


#### A couple of functions to define the infection rate

beta.fun = function(t,p,more){
    return( p['b0'] + p['b1']*sin(2*pi*t) + p['b2']*cos(2*pi*t) )
}

beta.dfdp = function(t,p,more){
    dfdp =  cbind(rep(1,length(t)), sin(2*pi*t), cos(2*pi*t)) 
    colnames(dfdp) = c('b0','b1','b2')
    return(dfdp)
}

betamore = list(beta.fun=beta.fun,
                beta.dfdp=beta.dfdp,
                beta.ind=c('b0','b1','b2'))


# Create a collocation basis to represent the state vector

rr = range(SEIRtimes)
knots = seq(rr[1],rr[2],2/52)
norder = 3
nbasis = length(knots)+norder-2

bbasis = create.bspline.basis(range=rr,norder=norder,nbasis=nbasis,breaks=knots)

# To get an initial estimate of the states we smooth the observed I component
# and set the other coefficients to zero.  

DEfd = smooth.basis(SEIRtimes,logdata[,3],fdPar(bbasis,1,0.1))

plotfit.fd(log(SEIRdata),SEIRtimes,DEfd$fd)

coefs = cbind(matrix(0,bbasis$nbasis,2),DEfd$fd$coefs)
DEfd = fd(coefs,bbasis)


# We will want to represent the state variables on the log scale so that they remain
# always positive. To do this, we set the 'posproc' component of LS.setup to 1. We
# will also compare the logstate to the log of the data directly, and therefore set
# 'poslik' to 0.

# We call LS.setup first and use the outputted lik and proc objects to pull the 
# coefficients for the unobserved state variables into line with the 
# differential equation. 

objs = LS.setup(SEIRpars,fn=SEIRfn,fd.obj=DEfd,more=betamore,data=data,times=SEIRtimes,
  posproc=TRUE,poslik=FALSE,names=SEIRvarnames,lambda=c(100,1,1))

proc = objs$proc
lik = objs$lik

res1 = FitMatchOpt(coefs=coefs,which=1:2,proc=proc,pars=SEIRpars,meth='nlminb')


# Let's have a look at the result

DEfd1 = fd(res1$coefs,bbasis)
plot(DEfd1,ylim=c(5,13))
points(SEIRtimes,logdata[,3])

# We can now run an initial smooth using the estimated coefficients as starting points. 

res2 = inneropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=proc,lik=lik,coefs=res1$coefs,in.meth='nlminb')

# Has this changed much?

DEfd2 = fd(res2$coefs,bbasis)

plot(DEfd2,lwd=2,ylim=c(5,13))
lines(DEfd1)


# And we can call the optimizing functions. 

res3 = outeropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=proc,lik=lik,coefs=res2$coefs,
  active=c('i','b0','b1','b2'))


# Some plots

## First, the estimated trajectories

DEfd3 = fd(res3$coefs,bbasis)
plot(DEfd3,lwd=2,ylim=c(5,14))


# Let's compare this to the data

plotfit.fd(logdata[,3],SEIRtimes,DEfd3[3],ylab='Fit to Data')


# We can also look at the discrepancy between the estimated trajectory and the
# differential equation

traj = eval.fd(SEIRtimes,DEfd3)     # estimated trajectory
colnames(traj) = SEIRvarnames

dtraj = eval.fd(SEIRtimes,DEfd3,1)  # derivative of the estimated trajectory

ftraj = proc$more$fn(SEIRtimes,traj,res3$pars,proc$more$more)   # Trajectory predicted by ODE


X11()
matplot(SEIRtimes,dtraj,type='l',lty=2,ylim =c(-10,10),ylab='SEIR derivatives' )
matplot(SEIRtimes,ftraj,type='l',lty=1,add=TRUE)

X11()
matplot(SEIRtimes,dtraj-ftraj,type='l',ylim=c(-4,4),ylab='Fit to Model')





## The alternative is to exponentiate the state before we compare to original data.
# This can take a very long time and is only recommended if you really need to do
# it, or have a couple of hours to wait. 


objs2 = LS.setup(SEIRpars,fn=SEIRfn,fd.obj=DEfd,more=betamore,data=data,times=SEIRtimes,
  posproc=TRUE,poslik=TRUE,names=SEIRvarnames,SEIRparnames,lambda=c(100,1,1))

lik2 = objs2$lik
proc2 = objs2$proc

res2 = inneropt(data=data,times=SEIRtimes,pars=res3$pars,proc=proc2,lik=lik2,coefs=res3$coefs)

res3 = outeropt(data=data,times=SEIRtimes,pars=res3$pars,proc=proc2,lik=lik2,coefs=res3$coefs,
  active=c('i','b0','b1','b2'))



###############################################################################
# Some more basic setup operations
#
# Here we go through the steps necessary to set up the SEIR equations manually. 
# This allows us several options for dealing with the positivity of the 
# state vector. 
#
# 1. We can ignore it and hope for the best. 
#
# 2. We can take a log transform of the ODE and then exponentiate the solutions
# to compare to the data
#
# 3. We can take a log transform of the ODE and compare this to the log data. 
###############################################################################


# First of all, we need the values of the basis at the observation times and the 
# quadrature points. 


qpts = 0.5*(knots[1:(length(knots)-1)] + knots[2:length(knots)])

bvals.obs = Matrix(eval.basis(times,bbasis),sparse=TRUE)

bvals = list(bvals = Matrix(eval.basis(qpts,bbasis),sparse=TRUE),
            dbvals = Matrix(eval.basis(qpts,bbasis,1),sparse=TRUE))


### This proc object is just the standard squared error deviation
# from the right hand side of a differential equation

sproc = make.SSEproc()
sproc$bvals = bvals
sproc$more = make.SEIR()
sproc$more$more = betamore
sproc$more$qpts = qpts
sproc$more$weights = matrix(1,length(qpts),3)%*%diag(c(1e2,1e0,1e0))
sproc$more$names = SEIRvarnames
sproc$more$parnames = SEIRparnames


### However, ODEs are often much more numerically stable if represented on a
# log scale. The make.logtrans() function will take the right hand side functions
# and derivatives defined for any differential equation and provide the equivalent
# system for the log state. Note that this does affect the way  you represent
# your state when considering the observations. 

lsproc = make.SSEproc()
lsproc$bvals = bvals
lsproc$more = make.logtrans()
lsproc$more$more = make.SEIR()
lsproc$more$more$more = betamore
lsproc$more$qpts = qpts
lsproc$more$weights = matrix(1,length(qpts),3)%*%diag(c(1e2,1e0,1e0))
lsproc$more$names = SEIRvarnames
lsproc$more$parnames = SEIRparnames

### Lik objects, this is the standard squared error. 

slik = make.SSElik()
slik$bvals = eval.basis(times,bbasis)
slik$more = make.id()
slik$more$weights = array(1,dim(data))
slik$more$names = SEIRvarnames
slik$more$parnames = SEIRparnames

# Log transform transformation for the lik object. For this, we note that we 
# have represented the trajectory on the log scale and will need to transform 
# back, this makes the numerics much harder and it can take a very long time to 
# converge.  

lslik = make.logstate.lik()
lslik$bvals = slik$bvals
lslik$more$weights = slik$more$weights
lslik$more = slik
lslik$more$parnames = SEIRparnames


# Numerically things work much better on the log scale

dfd = data2fd(logdata[,3],times,bbasis)

coefs = matrix(0,nbasis,3)
coefs[,3] = dfd$coefs
 
res = FitMatchOpt(coefs=coefs,which=1:2,proc=lsproc,pars=SEIRpars,meth='nlminb')

res2 = inneropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=lsproc,lik=lslik,coefs=res$coefs)

res3 = outeropt(data=data,times=SEIRtimes,pars=SEIRpars,proc=lsproc,lik=lslik,coefs=res2$coefs,
  active=c('i','b0','b1','b2'))


## Or just log the observations, this is faster. 

res2 = inneropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=lsproc,lik=slik,coefs=res$coefs)

res3 = outeropt(data=logdata,times=SEIRtimes,pars=SEIRpars,proc=lsproc,lik=slik,coefs=res2$coefs,
  active=c('i','b0','b1','b2'))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Groundwater level measurements as a function of rainfall.  
#
# The data were collected in North Vancouver near a site of a 
# mud slide that caused two fatalities in 2004.
#
# The groundwater levels are a relatively smooth process with a
# rather high signal to noise ratio.  Rainfall events, however,
# are sporadic and can be considered as events with times and 
# intensities, called a marked point process in statistical theory.
#
# Groundwater level is modeled by a first order linear differential
# equation forced by rainfall.
#
# A constant coefficient version of the equation is
# 
# $$\dot{G}(t) = - \beta G(t) + \alpha R(t - \delta) + \mu$$
#
# and a variable coefficient version is
#
# $$\dot{G}(t) = - \beta(t) G(t) + \alpha(t) R(t - \delta) + \mu(t)$$
#
# Coefficient \beta controls the speed of the response
# of groundwater to a rainfall event.  
#
# Coefficient \alpha controls the size of the response
# of groundwater to a rainfall event.   The gain in the system
# is defined as K = \alpha / \beta.
#
# Constant \mu is an adjustment to the overall level of
# groundwater required because the origin of the measurements
# can be viewed as arbitrary or irrelevant.
#
# Lag constant \delta allows a time lapse before groundwater
# begins to respond to a rainfall event.  It was estimated by
# inspecting plots to be three hours for the purposes of these
# analyses, and is therefore not actually estimated by these
# analyses themselves.
#
# The following code allows for each of the remaining three
# parameters to be either constant, or themselves functions 
# of time.  That is,
#
# Although the constant model does a fair job of 
# describing the data, allowing for a mild time-variation in
# these parameters substantially improved the fit.  We used
# five order 4 B-splines to model this time dependency.  
#
# Last modified 26 April 2010
#
## Setting up the data 
#

# The analyses requires access to the profiled estimation software.

library("CollocInfer")

# Input the data to be analyzed.

#load("groundwater")
#load("rainfall")

#  rainfall should be lagged by about three hours

#Ninput = length(rainfall)
#lag    = 3
#rainfall = c(rep(0,2*lag), rainfall[1:(Ninput-2*lag)])

#  a subset of the data are selected, running from day 79 into day 92
#  if all the data are used, R has to allocate additional storage

# index = 1597:3024  
#index = 1907:2221  #  here we select a smaller set of data 
#yobs = as.matrix(groundwater[index])
#zobs = as.matrix(rainfall[index])
#N = length(yobs)
#tobs = (index - index[1])


yobs = NSgroundwater
zobs = NSrainfall
tobs = NStimes

N    = length(tobs)
rangeval = c(0,N-1)



#  plot the groundwater data

par(mfrow = c(2,1))
plot(tobs, yobs, "p",  xlab="Time (hours)", 
                       ylab="Groundwater level (metres)")
plot(tobs, zobs, "p",  xlab="Time (hours)", 
                       ylab="Hourly rainfall (millimetres)")

#  set up cumulative rainfall

#  Set up a basis for rainfall itself: 
#     An order 1 or piece-wise constant spline

norder = 1
nbasis = N + norder - 2
rainbasis  = create.bspline.basis(rangeval, nbasis, norder)
rainfd = smooth.basis(tobs, zobs, rainbasis)$fd

#  plot the data

plotfit.fd(zobs, tobs, rainfd, nfine=N)

## The basis for representing groundwater: G(t)

# Define the knot sequence.  A knot is placed at each mid-hour
# point, giving maximum capacity to fit the data.  Variable
# rangeval is a vector of length 2 containing the lower and
# upper time limits of the observations, and is in file 
# NorthShoreData.mat


knots  = c(rangeval[1],
           seq(rangeval[1]+0.5, rangeval[2]-0.5, len=N-1), 
           rangeval[2])
norder   = 3
nbasis   = length(knots) + norder - 2
basisobj = create.bspline.basis(rangeval, nbasis, norder, knots)
 
# The constant coefficient basis is required for constant beta and alpha.  

conbasis = create.constant.basis(rangeval)

## Definition of functional basis objects for the coefficient functions

# Set up the basis for the coefficient function \beta(t) or
# \beta(G) determining the speed of response and recovery
#
# select the command that specifies the argument of \beta
#
# This basis is for \beta(t) as a function of time.  You may
# change the number of basis functions that are used.
#
# Here we use a constant basis for each coefficient, but at end
# we can do that analysis again with five order 5 B-spline basis functions.

nbetabasis = 1
#nbetabasis = 5
if (nbetabasis > 1) {
    betabasis = create.bspline.basis(rangeval, nbetabasis)
} else {
    betabasis = conbasis
}

# Set up the basis for the coefficient function \alpha(t)
# multiplying rainfall

nalphabasis = 1
#nalphabasis = 5
if (nalphabasis > 1) {
    alphabasis = create.bspline.basis(rangeval, nalphabasis)
} else {
    alphabasis = conbasis
}

#  store the bases and rainfall object in a list "more"

more = vector("list",0)
more$betabasis  = betabasis
more$alphabasis = alphabasis
more$rainfd     = rainfd

# Definition of the required named list containing the right side
#  evaluation functions.  

NSfn = make.NS()

## Prepare the constant coefficient analysis

# Set up the functional parameter object that defines the
# roughness penalty.  Here an important choice must be made:
# how large should the smoothing parameter \lambda be?
# One may be advised to begin with a lowish value, such as
# \lambda = 1, complete the analysis with this value, then
# increase the value to something like \lambda = 100, using
# as starting values for the parameter array pars at this stage
# the array newpars computed at the profiled estimation analysis
# for the previous smaller value.  

lambdaDE = 1e0  # A good value for the initial analysis.     
penorder = 1    # The penalty is first order
GfdPar   = fdPar(basisobj,penorder,lambdaDE)     

# Smooth the data with this roughness penalty to get initial
# values for coefficients defining G(t).

DEfd0 = smooth.basis(tobs, yobs, GfdPar)$fd

# plot the data and the fit

plotfit.fd(yobs, tobs, DEfd0, 
           xlab="Time (hours)", ylab="Groundwater level (metres)", 
           title="")

# extract the coefficient array

coefs0 = DEfd0$coefs

# Initial parameter values for constant bases

pars0 = matrix(0, nbetabasis + nalphabasis, 1) 

#pars0 = matrix(c(rep(res1$pars[1],nbetabasis), 
#                 rep(res1$pars[1],nalphabasis)), 
#                 nbetabasis + nalphabasis, 1) 

## Profiled estimation step

control.out = list()
control.out$trace  = 6
control.out$tol    = 1e-6

# This is the main optimization for profiled estimation.  
# It uses initial values for parameters and coefficients that
# are set up above.  This requires a lot of computation, and
# will take a number of minutes per iteration.  Find something
# else to do in the meantime.

lambda = 1e0

#  estimate the parameters using differencing

res0 = Profile.LS(make.NS()$fn, yobs, tobs, pars0, coefs0, basisobj, 
                  lambda, more = more, out.meth='ProfileGN', 
                  control.out=control.out)

res0$pars  #  display the parameter estimates

res0$outer.result$value  # display the final function value


#  estimate the parameters using the derivative evaluation functions

lambda = 1e0

res1 = Profile.LS(NSfn, yobs, tobs, pars0, coefs0, basisobj, lambda,
                  more=more, out.meth='ProfileGN', 
                  control.out=control.out)

res1$pars  #  display the parameter estimates

res1$outer.result$value  # display the final function value

#  set up the functional data object for the fit

DEfd1 = fd(res1$coefs, basisobj)

#  increase lambda a few times, starting each time with 
#  parameters estimated from the last value

lambda = 1e2

res2 = Profile.LS(NSfn, yobs, tobs, res1$pars, res1$coefs, basisobj, lambda,
                  more = more, out.meth='ProfileGN', 
                  control.out=control.out)

res2$pars  #  display the parameter estimates

DEfd2 = fd(res2$coefs, basisobj)

lambda = 1e4

res3 = Profile.LS(NSfn, yobs, tobs, res2$pars, res2$coefs, basisobj, lambda,
                  more = more, out.meth='ProfileGN', 
                  control.out=control.out)

res3$pars  #  display the parameter estimates

DEfd3 = fd(res3$coefs, basisobj)

lambda = 1e6

res4 = Profile.LS(NSfn, yobs, tobs, res3$pars, res3$coefs, basisobj, lambda,
                  more = more, out.meth='ProfileGN', 
                  control.out=control.out)

res4$pars  #  display the parameter estimates

DEfd4 = fd(res4$coefs, basisobj)

lambda = 1e8

res5 = Profile.LS(NSfn, yobs, tobs, res4$pars, res4$coefs, basisobj, lambda,
                  more = more, out.meth='ProfileGN', 
                  control.out=control.out)

res5$pars  #  display the parameter estimates

DEfd5 = fd(res5$coefs, basisobj)

#  compare fits to the data

par(mfrow=c(1,1))
plotfit.fd(yobs, tobs, DEfd1, 
           xlab="Time (hours)", ylab="Groundwater level (metres)", 
           title="")
lines(DEfd2)
lines(DEfd3)
lines(DEfd4)
lines(DEfd5, lty=2)

#  By the time we get to res4, the fit to the data is essentially
#  a soluton to the differential equation.  We see that it has
#  rather less capacity to fit the data, but still captures some
#  of the shapes in the data.

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load some things in 

library('CollocInfer')

# The Chemostat equations represent a four-species Chemostat plus the resource 
# of Nitrogen. There are two species of Algae with varying defenses against
# Rotifers. The Rotifers themselves are divided into two class -- breeding
# and senescent, although these two are very tightly coupled. 
#
# A full description of these equations can be found in the user manual. 
# The five state variables for the equations are
# 
# N - nitrogen content in the Chemostat
# C1 - Algal type 1
# C2 - Algal type 2 
# B - Breeding Rotifers
# S - Senescent Rotifers
#
# The system has 16 parameters, also described in the user manual. Notable
# features include that only the sums C1+C2 and B+S can be observed. Further, 
# an unknown fraction of each is counted at each time. This requires us to 
# set up a model for the observation process along with the ODE. 


# First we load up some data

data(ChemoData)

# The first two of these parameters give the fractions of Algae and Rotifers 
# that are counted. The remaining parameters are all positive and using 
# their logged values is helpful. 

logpars=c(ChemoPars[1:2],log(ChemoPars[3:16]))

# Parameters 'p1' and 'p2' represent relative palatability of the two algal
# clones, as such only one can be estimated and we fix p2 = 0. 

active = c(1:2,5,7:16)     

# We'll choose a fairly large value of lambda. 

lambda = rep(200,5)  

# We need some basis functions

rr = range(ChemoTime)
knots = seq(rr[1],rr[2],by=0.5)
bbasis = create.bspline.basis(rr,norder=4,breaks=knots)

# We will also have to set up the basis matrices manually. 

mids = c(min(knots),(knots[1:(length(knots)-1)] + 0.25),max(knots))

bvals.obs = eval.basis(ChemoTime,bbasis)

bvals.proc = list(bvals = eval.basis(mids,bbasis),
                  dbvals = eval.basis(mids,bbasis,1));


# We can now set up the proc object. We will want to take a log transformation
# of the state here for numerical stability. In general it is better to do 
# finite differencing AFTER the log transformation rather than before it. 

proc = make.SSEproc()                    # Sum of squared errors
proc$bvals = bvals.proc                  # Basis values
proc$more = make.findif.ode()            # Finite differencing
proc$more$more = list(fn=make.logtrans()$fn,eps=1e-8) # Log transform
proc$more$more$more = list(fn=chemo.fun) # ODE function
proc$more$qpts = mids                    # Quadrature points
proc$more$weights = rep(1,5)*lambda      # Quadrature weights (including lambda)
proc$more$names = ChemoVarnames          # Variable names
proc$more$parnames = ChemoParnames       # Parameter names


# For the lik object we need to both represent the linear combination transform
# and we need to model the observation process. 

# First to represent the observation process, we can use the genlin
# functions. These produce a linear combination of the the states 
# (they can be used in proc objects for linear systems, too).  

temp.lik = make.SSElik()
temp.lik$more = make.genlin()

# Genlin requires a more object with two elements. The 'mat' element
# gives a template for the matrix defining the linear combination. This is
# all zeros 2x5 in our case for the two observations from five states. 
# The 'sub' element specifies which elements of the parameters should be
# substituted into the mat element. 'sub' should be a kx3 matrix, each
# row defines the row (1) and column (2) of 'mat' to use and the element
# of the parameter vector (3) to add to it. 

temp.lik$more$more = list(mat=matrix(0,2,5,byrow=TRUE), 
  sub = matrix(c(1,2,1,1,3,1,2,4,2,2,5,2),4,3,byrow=TRUE))
temp.lik$more$weights = c(100,1)

# Finally, we tell CollocInfer that the trajectories are represented on
# the log scale and must be exponentiated before comparing them to the data. 

lik = make.logstate.lik()
lik$more = temp.lik
lik$bvals = bvals.obs


# Now lets try running this

# Because we don't have direct observations of any state, we'll use a starting 
# smooth obtained from generating some ODE solutions

y0 = log(c(2,0.1,0.4,0.2,0.1))
names(y0) = ChemoVarnames

odetraj = lsoda(y0,ChemoTime,func=chemo.ode,logpars)

DEfd = smooth.basis(ChemoTime,odetraj[,2:6],fdPar(bbasis,int2Lfd(2),1e-6))
C = DEfd$fd$coef


# Now, with parameters fixed, we'll estiamte coefficients. 

res = inneropt(coefs=C,pars=logpars,times=ChemoTime,data=ChemoData,
  lik=lik,proc=proc,in.meth='optim')

# We'll for the trajectory and also the appropriate sum of exponentiated
# states to compare to the data. 

C = matrix(res$coefs,dim(C))
traj = lik$bvals%*%C
obstraj = lik$more$more$fn(ChemoTime,exp(traj),logpars,lik$more$more$more)

# Plot these against the data

X11()
par(mfrow=c(2,1))
plot(obstraj[,1],type='l',ylab='Chlamy',xlab='',cex.lab=1.5,cex.axis=1.5)
points(ChemoData[,1])
plot(obstraj[,2],type='l',ylab='Brachionus',xlab='days',cex.lab=1.5,cex.axis=1.5)
points(ChemoData[,2])


# Now we can continue with the outer optimization

res2 = outeropt(pars=logpars,times=ChemoTime,data=ChemoData,coef=C,
    lik=lik,proc=proc,active=active,in.meth='optim',out.meth='nlminb')

# We'll extract the resulting parameters and coefficients. 

npars = res2$pars
C = as.matrix(res2$coefs,dim(C))

# And obtain an estimated trajectory and the exponentiated sum to comprare
# to the data. 

traj = lik$bvals%*%C
ptraj = lik$more$more$fn(ChemoTime,exp(traj),npars,lik$more$more$more)

# Lets have a look at how much we changed our parameters on the original 
# scale. 

new.pars = npars
new.pars[3:16] = exp(new.pars[3:16])

print(ChemoPars)
print(new.pars)
print(new.pars/ChemoPars)

# Now we can produce a set of diagnostic plots. 

# Firstly, a representation of the trajectory compared to the data. 

X11()
par(mfrow=c(2,1))
plot(ChemoTime,ptraj[,1],type='l',ylab='Chlamy',xlab='',cex.lab=1.5)
points(ChemoTime,ChemoData[,1])
plot(ChemoTime,ptraj[,2],type='l',ylab='Brachionus',xlab='days',cex.lab=1.5)
points(ChemoTime,ChemoData[,2])


# Now we'll plot both the derivative of the trajectory and the value of the
# differential equation right hand side at each point. This represents the 
# fit to the model. 

traj2 = proc$bvals$bvals%*%C
dtraj2 = proc$bvals$dbvals%*%C

colnames(traj2) = ChemoVarnames
ftraj2 = proc$more$fn(proc2$more$qpts,traj2,npars,proc$more$more)


X11()
par(mfrow=c(5,1),mai=c(0.3,0.6,0.1,0.1))
for(i in 1:5){
  plot(mids,dtraj2[,i],type='l',xlab='',ylab=ChemoVarnames[i],
    cex.lab=1.5,ylim=c(-0.5,0.5))
  lines(mids,ftraj2[,i],col=2,lty=2)
  abline(h=0)
}
legend('topleft',legend=c('Smooth','Model'),lty=1:2,col=1:2,cex=1.2)

# Solving the differential equation from the estiamted initial condition
# of the trajectory allows us to compare the qualitative behavior of 
# our estimate to that of the differential equation. 

y0 = traj[1,]
names(y0) = ChemoVarnames
odetraj = lsoda(y0,ChemoTime,func=chemo.ode,parms=npars)


X11()
par(mfrow=c(2,1))
matplot(ChemoTime,traj,col=1,type='l',lwd=3,cex.lab=1.5,cex.axis=1.5,
    ylab='',cex.main=1.5,main='Reconstructed Trajectories')
legend(x='topright',legend=ChemoVarnames,lwd=3,col=1,lty=1:5)
matplot(ChemoTime,odetraj[,2:6],col=1,type='l',lwd=3,cex.axis=1.5,cex.lab=1.5,
    ylab='',cex.main=1.5,main='ODE Solution')

# We can also compare the pattern of observations predicted by the differential
# equation and that estimated by our methods. 

otraj = lik$more$more$fn(ChemoTime,exp(odetraj[,2:6]),npars,lik$more$more$more)

X11()
par(mfrow=c(2,1))
matplot(ChemoTime,ptraj,type='l',lwd=2,xlab='days',cex.lab=1.5,ylab='',
  cex.axis=1.5,cex.main=1.5,main='Predicted Observations -- Smooth')
matplot(ChemoTime,ChemoData,add=TRUE,pch = c(1,2))
legend('topright',legend=c('Algae','Rotifers'),pch=1:2,col=1:2)
matplot(ChemoTime,otraj,type='l',lwd=2,xlab='days',cex.lab=1.5,ylab='',
  cex.axis=1.5,cex.main=1.5,main='Predicted Observations -- ODE')
legend('topright',legend=c('Algae','Rotifers'),lty=1:2,col=1:2,lwd=2)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

Ires1	= Smooth.LS(make.Henon(),data=Y,times=t,pars=hpars2,coefs,basisvals=NULL,
  lambda=lambda,in.meth='nlminb',control.in=control.in,pos=FALSE,discrete=TRUE)
  
Ores1 = Profile.LS(make.Henon(),data=Y,t,pars=hpars2,coefs,basisvals=NULL,
  lambda=lambda,in.meth='nlminb',out.meth='nls',control.in=control.in,
  control.out=control.out,pos=FALSE,discrete=TRUE)

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
  var=var,in.meth='nlminb',control.in=control.in,pos=FALSE,discrete=TRUE)

Ores3 = Profile.multinorm(fn=make.Henon(),data=Y,times=t,pars=hpars2,coefs=coefs,
  basisvals=NULL,var=var,fd.obj=NULL,more=NULL,quadrature=NULL,in.meth='nlminb',
  out.meth='BFGS',control.in=control.in,control.out=control.out,eps=1e-6,
  active=NULL,pos=FALSE,discrete=TRUE)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Example Diagnostics -- Learning the FitzHugh-Nagumo Equations
#
# Demonstration estimation of function functions. This also provides
# a template for random-effects treatment within profiling.

library(fda)
library(odesolve)
library(MASS)
library('CollocInfer')

source('../Devel/Diagnostics.R')

# First Create Some Data

t = seq(0,20,0.05)

pars = c(0.2,0.2,3)
names(pars) = c('a','b','c')

x0 = c(-1,1)
names(x0)= c('V','R')
y = lsoda(x0,t,make.fhn()$fn.ode,pars)
y = y[,2:3]

data = y + 0.2*matrix(rnorm(802),401,2)

# Now a basis object

knots = seq(0,20,0.2)
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	norder=norder,breaks=knots)


# Initial values for coefficients will be obtained by smoothing

DEfd = data2fd(data,t,bbasis)
coefs = DEfd$coefs

names(coefs) = c('V','R')

# Usual meta-parameters; quadrature points, weights and knots

lambda = c(100,100)
qpts = knots
qwts = rep(1/length(knots),length(knots))

qwts = qwts%*%t(lambda)
weights = array(1,dim(data))


# Now I define a measurement process log likelihood along with some
# additional features: in this case it's squared error. 

varnames = c('V','R')
parnames = c('a','b','c')


likmore = make.id()
likmore$weights = weights

lik = make.SSElik()
lik$more = likmore
lik$bvals = eval.basis(t,bbasis)


# Proc is a process log likelihood -- in this case treated as squared
# discrepancy from the ODE definition. 

procmore = make.genlin()
procmore$names = varnames
procmore$parnames = parnames
procmore$more = list(mat=matrix(0,2,2),sub= matrix(c(1,1,1,1,2,2,2,1,3,2,2,4),4,3,byrow=TRUE))


procmore$weights = qwts
procmore$qpts = qpts


proc = make.SSEproc()
proc$more = procmore
proc$bvals = list(bvals=eval.basis(procmore$qpts,bbasis,0),
		dbvals = eval.basis(procmore$qpts,bbasis,1))



spars = c(0,1,-1,0)

Ires = inneropt(data,times=t,spars,coefs,lik,proc,in.meth='nlminb')

Ores = outeropt(data=data,times=t,pars=spars,coefs=Ires$coefs,lik=lik,proc=proc,
  in.meth="nlminb",out.meth="nlminb")



traj = as.matrix(proc$bvals$bvals %*% Ores$coefs)
dtraj = as.matrix(proc$bvals$dbvals %*% Ores$coefs)
ftraj = dtraj - proc$more$fn(proc$more$qpts,dtraj,Ores$pars,proc$more$more)



par(mfrow=c(2,2))
for(i in 1:2){
  for(j in 1:2){
    plot(traj[,i],ftraj[,j],type='l')
  }
}




## Now we estimate some forcing functions

fbasis = create.bspline.basis(range=range,nbasis=23,norder=4)


dproc = make.SSEproc()
dproc$more = make.diagnostics()
dproc$more$qpts = procmore$qpts
dproc$more$weights = procmore$weights
dproc$more$more = procmore
dproc$more$more$p = Ores$pars
dproc$more$more$which = 1:2
dproc$more$more$psi = eval.basis(procmore$qpts,fbasis)
dproc$bvals = list(bvals=eval.basis(procmore$qpts,bbasis,0),
		dbvals = eval.basis(procmore$qpts,bbasis,1))


dpars = rep(0,2*fbasis$nbasis)  
  
dOres = outeropt(data=data,times=t,pars=dpars,coefs=Ires$coefs,lik=lik,proc=dproc,
  in.meth="nlminb",out.meth="nlminb")  


# Trajectories

force = dproc$more$more$psi %*% matrix(dOres$par,fbasis$nbasis,2)
traj = dproc$bvals$bvals %*% dOres$coefs

plot(traj[,1],force[,1],type='l')



