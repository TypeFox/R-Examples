##########
# This demonstration illustrates the use of log transformations with
# CollocInfer. 
#  We use the SEIR equations, various modifications of which are used to model
#  the spread of disease in epidemeiology.  The letters stand for:
#
#  Susceptibles:  the large subpopulation of those who have not been infected
#  Exposed:       the smaller subpopulation of those who have been infected,
#                 but are not yet infectious
#  Infectious:    those who have already been infected, and are now also 
#                 able to infect others
#  Recovered:     these were infected, but are now recovered and are, usually,
#                 immune to further infections.
#
#  The equations are used to model measles infections in the Canadian province
#  of Ontario over a two-year period.  The data analyzed here are derived from
#  a larger dataset analyzed in greater detail by Hooker, Ellner, De Vargas
#  Roditi and Earn (2011).  
#
#  The basic equations are here modified to all for a seasonally varying
#  infection rate, a realistic feature for modeling childhood diseases where
#  the infection rate increases during the school year.  The function that
#  provides the seasonal variability is beta(t).
#
# The equations are as follows.  Each equation has been organized so that
# the first term on the right defines the exponential decline of the 
# subpopulation if it were left to itself, and the remaining terms, called
# "forcing" terms, define the inputs from the other variables and external
# sources that cause the population to increase (or decline further).
#
# Sdot = -nu*S         + mu - beta(t)*S*(I+i) (Susceptibles)
# Edot = -(sigma+nu)*E + beta(t)*S*(I+i)      (Exposed)
# Idot = -(gamma+nu)*I + sigma*E              (Infectious)
# Rdot = -nu*R         + gamma*I               (Recovered)
# 
# Here beta(t) is parameterized by a sinusoidal function plus a constant,
#      beta(t) = p_1 + p_2 sin(2 pi t) + p_3 cos(2 pi t)
# where the period is one year.
#
# Other parameters in the equations are
# mu    - the arrival rate of new individuals into the population at large
# nu    - the death rate in the population at large
# i     - an external source of infected individuals
# sigma - the rate of movement from Exposed to Infectious.
# gamma - the rate of recovery from infection.
#
# Because the number of susceptibles is usually many times the number of
# infected individuals, S is on a quite different scale than E and I. 
# In these data, S varies around about 450,000, but E and I vary between
# about 400 and 9000. 
# Moreover, all of these variables cannot in preinciple be negative.  Both 
# computational efficiency and positivity argue for fitting the logarithm of the 
# data, and also modelling the logarithm of each variable. 
# This requires a simple modification of the SEIR equations above, consisting
# of dividing the right side of the equations by the exponential of the 
# respective variables.
# The CollocInfer package contains a number of useful tools that let you make 
# this transition without needing to re-code all of your differential equations.


library('CollocInfer')
                        
# The following  objects are a part of the package:

SEIRvarnames # the names of the variables
SEIRparnames # the names of the parameters

SEIRtimes    # the times of observation
SEIRdata     # the number of infected individuals
SEIRpars     # the initial estimates of parameter values

# Now augment the data for I  so that S and E measurements are listed as NA

data = cbind(matrix(NA,length(SEIRdata),2),SEIRdata)

# Calculate the natural logs of the observations  
# Note: using common logs would have made the results easier for clients to
# interpret, and would have required dividing the SEIR right sides above by
# 10 to the power of the variable values rather than vy their exponentials.

logdata = log(data)

# Function make.SEIR is distributed with the package, and defines 
# the right sides of the differential equations as well as various of their
# partial derivatives required during the computation.   The code can be 
# viewed by typing "make.SEIR" into the command window in R.

SEIRfn = make.SEIR()

# Define the time-varying infection rate beta  as a sinusoid with period
# one year plus a constant

beta.fun = function(t,p,more){
    return( p['b0'] + p['b1']*sin(2*pi*t) + p['b2']*cos(2*pi*t) )
}

# Define as well its partial derivatives with respect to the coefficients

beta.dfdp = function(t,p,more){
    dfdp =  cbind(rep(1,length(t)), sin(2*pi*t), cos(2*pi*t)) 
    colnames(dfdp) = c('b0','b1','b2')
    return(dfdp)
}

#  Bundle these two functions as well as the names of their coefficients
#  into a list object called betamore

betamore = list(beta.fun =beta.fun,
                beta.dfdp=beta.dfdp,
                beta.ind =c('b0','b1','b2'))

# Set up a B-spline functional basis object to represent the state vector

rr     = range(SEIRtimes)       #  the range of observations times
knots  = seq(rr[1],rr[2],2/52)  #  knots at 52 equally spaced values
norder = 3                      #  the order of the B-spline basis functions,
                                #  in this case piece-wise quadratic
nbasis = length(knots)+norder-2 #  the number of basis functions

#  set up the basis object

bbasis = create.bspline.basis(range=rr, norder=norder, nbasis=nbasis, 
                              breaks=knots)

# To get an initial estimate of the states we smooth the observed I component
# and set the other coefficients to zero.  

# smooth the log I values

DEfd = smooth.basis(SEIRtimes,logdata[,3],fdPar(bbasis,1,0.1))

# plot the smooth plus data

plotfit.fd(log(SEIRdata),SEIRtimes,DEfd$fd)

# Augment the coefficient vector by two columns.  The first column
# contains all 13's, which is about the logarithm of the susceptible 
# population size, and the second contains all 7's, about the estimated
# size of the exposed subpopulation.  The ratio of the two population sizes
# is here approximated by exp(13-7), or about 400.

coefs = cbind(matrix(13,bbasis$nbasis,1),
              matrix( 7,bbasis$nbasis,1),
              DEfd$fd$coefs)

# set up the functional data object for the three variables

DEfd = fd(coefs,bbasis)

# We use the usual least squares measures of the fidelity of the estimated
# functions S(t), E(t) and I(t) to the data and to the differential equation.
# This means that we can save the user the effort of setting up the lik and
# proc objects by calling function LS.setup, which output these objects.

# We will want to represent the state variables on the log scale so that they 
# remain always positive. To do this, we set the 'posproc' argument of 
# LS.setup to 1. We will also compare the log I(t) to the log of the data 
# directly, and therefore set argument 'poslik' to 0, the default.

#  Here, too, we specify the level of emphasis on fitting each differential
#  equation by setting the corresponding smoothing or bandwidth parameter
#  lambda.  The following two statements specify the lambda to be 100 for
#  S, and 1 for E and I.

SEIRlambda = rep(1,3)
SEIRlambda[1] = 100

#  list object betamore is now passed into LS.setup, too, in order to make
#  it available as a functional parameter defined by its three coefficients

#  run LS.setup

objs = LS.setup(SEIRpars, fn=SEIRfn, fd.obj=DEfd, more=betamore, data=data, 
                times=SEIRtimes, posproc=1, poslik=0, names=SEIRvarnames, 
                lambda=SEIRlambda)
                
#  Set up the proc and lik objects for this problem
                
SEIRproc = objs$proc
SEIRlik  = objs$lik

#  The next task is to improve the initial coefficient estimates 13 and 7
#  defined above for S and E, respectively by using function FitMatchOpt
#  that optimize coefficients for unobserved variables.  The 'which' argument
#  of this function indicates which variables's coefficients are to be
#  optimized with respect to the fit to the observed data for I.
#  Function proc.time is used here to show the elapsed time, the time
#  titled 'user' is the current elapsed session time in seconds.
#  The optimization function to be used is function nlminb.


res1 = FitMatchOpt(coefs=coefs, which=1:2, proc=SEIRproc, pars=SEIRpars,
                   meth='nlminb')


# Let's have a look at the three functions and the fit to the I data

DEfd1 = fd(res1$coefs,bbasis)

plot(DEfd1,lwd=2,ylim=c(5,14),col=c(4,2,3),lty=c(1:3))
points(SEIRtimes,logdata[,3])  #  plot the data for I
legend(0,12,SEIRvarnames,col=c(4,2,3),lty=c(1:3),lwd=2)

# We can now run an initial smooth using the estimated coefficients as 
# starting points. 
# The optimization function to be used is function nlminb.

res2 = inneropt(data=logdata, times=SEIRtimes, pars=SEIRpars,
                proc=SEIRproc,lik=SEIRlik, coefs=res1$coefs, in.meth='nlminb')

#  define the new function estimates

DEfd2 = fd(res2$coefs,bbasis)

# Has this changed much?  Plot the new function estimates along with the old

plot(DEfd2,lwd=2,ylim=c(5,14),col=c(4,2,3),lty=c(1:3))
lines(DEfd1)
legend(0,12,SEIRvarnames,col=c(4,2,3),lty=c(1:3),lwd=2)

# All these analyses have been using the initial parameter estimates in
# order to get starting values for the coefficients defining the three
# functions.  Now we are ready to optimize the fit to the data with
# respect to the parameters in SEIRpars.  Here, however, we only optimize
# the values of the 'i' parameter defining the external input of infectives
# and the three coefficients defining the time-varying infection rate beta. 
# These are selected using the argument 'active' of outer optimization
# function outeropt.

# The computation time involved in this step is substantial (around 10 minutes
# on the computer used while preparing these notes).   You may
# want to disable output buffering by clicking the popup menu item 
# 'Buffered output'that you see when you click on the 'Misc' item at the
# top of the command window.


SEIRactive = c('i','b0','b1','b2')

res3 = outeropt(data=logdata, times=SEIRtimes, pars=SEIRpars,
                proc=SEIRproc, lik=SEIRlik, coefs=res2$coefs,
                active=SEIRactive)


#  display the initial and estimated parameter values

SEIRparsEst = res3$pars

print(matrix(c(SEIRpars[SEIRactive],SEIRparsEst[SEIRactive]),4,2))

# Plot the estimated trajectories

DEfd3 = fd(res3$coefs,bbasis)

plot(DEfd3,lwd=2,ylim=c(5,14),col=c(4,2,3),lty=c(1:3))
legend(0,12,SEIRvarnames,col=c(4,2,3),lty=c(1:3),lwd=2)


# We can now employ CollocInferPlots to explore fit to the data

out3 = CollocInferPlots(pars=res3$pars,coefs=res3$coef,lik=SEIRlik,proc=SEIRproc,data=logdata,times=SEIRtimes)


# This function just automates the following calculations:

# Compare the smooth of the data to the data itself:

plotfit.fd(logdata[,3],SEIRtimes,DEfd3[3],ylab='Fit to Data')

# 2: We can also look at the discrepancy between the estimated trajectory and the
# differential equation in terms of its derivatives

traj = eval.fd(SEIRtimes,DEfd3)     # estimated trajectory
colnames(traj) = SEIRvarnames

dtraj = eval.fd(SEIRtimes,DEfd3,1)  # derivative of the estimated trajectory

# derivative as modeled by the ODE by ODE

ftraj = SEIRproc$more$fn(SEIRtimes,traj,res3$pars,SEIRproc$more$more)   

#  Plot the derivatives of the trajectories along with the derivatives
#  defined by the right side.

X11()
matplot(SEIRtimes,dtraj,type='l',lty=2,ylim =c(-10,10),col=c(4,2,3),
        ylab='SEIR derivatives' )
matplot(SEIRtimes,ftraj,type='l',lty=1,add=TRUE,col=c(4,2,3))

#  3. Plot the differences between the actual and predicted derivatives

X11()
matplot(SEIRtimes,dtraj-ftraj,type='l',ylim=c(-4,4), col=c(4,2,3),
        ylab='Residuals')
        
        

## The alternative is to exponentiate the state before we compare to the data.
# This can take a very long time and is only recommended if you really need 
# to do it, or have a couple of hours to wait. 

objs2 = LS.setup(SEIRpars, fn=SEIRfn, fd.obj=DEfd, more=betamore, data=data,
                 times=SEIRtimes, posproc=1, poslik=1, 
                 names=SEIRvarnames, SEIRparnames, lambda=SEIRlambda)

SEIRlik2  = objs2$lik
SEIRproc2 = objs2$proc

res2 = inneropt(data=data, times=SEIRtimes, pars=res3$pars,
                proc=SEIRproc2, lik=SEIRlik2, coefs=res3$coefs)

res3 = outeropt(data=data, times=SEIRtimes, pars=res3$pars,
                proc=SEIRproc2, lik=SEIRlik2, coefs=res3$coefs,
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

bvals.obs = Matrix(eval.basis(SEIRtimes,bbasis),sparse=TRUE)

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
slik$bvals = eval.basis(SEIRtimes,bbasis)
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

DEfd = smooth.basis(SEIRtimes,logdata,fdPar(bbasis,1,0.5)) 

coefs = matrix(0,nbasis,3)
coefs[,3] = DEfd$fd$coefs[,3] 
 
res = FitMatchOpt(coefs=coefs, which=1:2, proc=lsproc, pars=SEIRpars,
                  meth='nlminb')

# This is the non-log scale run; once again, this can take a few hours

res2 = inneropt(data=data, times=SEIRtimes, pars=SEIRpars,
                proc=lsproc, lik=slik, coefs=res$coefs)

res3 = outeropt(data=data, times=SEIRtimes, pars=SEIRpars, proc=lsproc,
                lik=slik, coefs=res2$coefs, active=c('i','b0','b1','b2'))
                

## Or just log the observations, this is faster. 

res2 = inneropt(data=logdata, times=SEIRtimes, pars=SEIRpars,
                proc=lsproc, lik=slik, coefs=res$coefs)

res3 = outeropt(data=logdata, times=SEIRtimes, pars=SEIRpars,
                proc=lsproc, lik=slik, coefs=res2$coefs,
                active=c('i','b0','b1','b2'))
                