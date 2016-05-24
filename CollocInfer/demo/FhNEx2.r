# This demonstration file is intended to accompany the CollocInfer manual. 

library(CollocInfer)


# To ensure reproducibility
set.seed((2004*2007)/2014)

# The first section of this demo file provides the code to generate simulation data
# from the FitzHugh-Nagumo model. We will delete it later and call a standard set of data
# with data(FhNdata). This is for the purposes of providing reproducible results; the data we
# will use were generated with exactly the commands below. 


# First we will define some variable and parameter names

FhNvarnames = c('V','R')
FhNparnames = c('a','b','c')

# and initial conditions and parameters

x0 = c(-1,1)
names(x0) = FhNvarnames

FhNpars = c(0.2,0.2,3)
names(FhNpars) = FhNparnames


# The following is a function specifying the FitzHugh-Nagumo model in a form that
# lsoda works with 

fhn.ode <- function(times,x,p)
{
 dx = x;
 dimnames(dx) = dimnames(x);
 dx['V'] = p['c']*(x['V'] - x['V']^3/3 + x['R'])
 dx['R'] = -(x['V'] -p['a'] + p['b']*x['R'])/p['c']

 return(list(dx))
}
 
# We need the times at which we will observe the system
 
FhNplottimes = seq(0,20,0.05)

# And now we can create solutions to the equations

out = lsoda(x0,times=FhNplottimes,fhn.ode,FhNpars)

# and plot what the solutions look like

par(mar=c(5,5,1,1))
matplot(out[,1],out[,2:3],type='l',xlab='time',ylab='(V,R)',lwd=3,cex.lab=2.5,cex.axis=2.5)
legend('bottomleft',c('V','R'),lwd=3,col=1:2,lty=1:2,cex=1.5)

par(mar=c(5,5,1,1))
plot(out[,2],out[,3],type='l',lwd=3,xlab='V',ylab='R',cex.lab=2.5,cex.axis=2.5)

# We now add some noise to the values of the curves at a reduced set of sampling times:

FhNtimes = seq(0,20,0.5)
FhNn = length(FhNtimes)
y = lsoda(x0,times=FhNtimes,fhn.ode,FhNpars)
FhNdata = out[,2:3] + 0.06*matrix(2*FhNn,FhNn,2)


## Let's undo all this work and call a standard data set generated from this
# code with which to run the profiling proceedures. 

rm(list=ls())
data(FhNdata)


# In order to run the profiling proceedures, we need to define some objects. 

# The following code will define a basis

range = c(0,20)
knots = seq(0,20,0.5)
norder = 4

FhNbasis = create.bspline.basis(range=range,norder=norder,breaks=knots)

# And this will smooth the data

DEfd0 = smooth.basis(FhNtimes,FhNdata,fdPar(FhNbasis,int2Lfd(2),1))$fd

#  and produce a plot of how well this agrees with the smooth
par(mfrow=c(2,1),mar=c(5,5,2,1))
plotfit.fd(FhNdata, FhNtimes, DEfd0, cex.axis=2.5, cex.lab=2.5, lwd=2, cex=1.5)


# We can now extract the coefficients, which we will also require

coefs0 = DEfd0$coef
colnames(coefs0) = FhNvarnames

# CollocInfer requires the right hand side function to be defined somewhat differently
# to lsoda. Here we allow a matrix of values as input 

fhn.fun <- function(times,x,p,more)
{
 dx = x;
 dx[,'V'] = p['c']*(x[,'V'] - x[,'V']^3/3 + x[,'R'])
 dx[,'R'] = -(x[,'V'] -p['a'] + p['b']*x[,'R'])/p['c']
 return(dx)
}           


# Now we can choose a trade-off parameter and set up the objects that the profiling
# functions will use. 

lambda = 1000

profile.obj = LS.setup(pars=FhNpars,fn=fhn.fun,lambda=lambda,times=FhNtimes,
                       coefs=coefs0,basisvals=FhNbasis)
                       
proc = profile.obj$proc
lik = profile.obj$lik
                       
## Gradient matching can be obtained thr ParsMatchOpt and produces useful initial
# parameter estimates

Pres0  = ParsMatchOpt(FhNpars,coefs0,proc)
pars1 = Pres0$pars

## Inner optimisation to smooth the data using the differential equation as a penalty

Ires1 = inneropt(FhNdata,times=FhNtimes,pars1,coefs0,lik,proc)
coefs1 = Ires1$coefs

# Alternatively, Smooth.LS avoids the need to run LS.setup, and it also returs
# the lik and proc objects. 

Ires1.2 = Smooth.LS(fhn.fun,FhNdata,FhNtimes,FhNpars,coefs0,FhNbasis,lambda)

# And outer optimization to choose parameters

Ores2 = outeropt(FhNdata,FhNtimes,pars1,coefs1,lik,proc)

# Here the relevant objects to extract are parameter estimates and the corresponding
# coefficients. 

pars2 = Ores2$pars
coefs2 = Ores2$coefs

# Alternatively, Profile.LS will do both setup and profiling

Ores2.2 = Profile.LS(fhn.fun,FhNdata,FhNtimes,FhNpars,coefs0,FhNbasis,lambda)


# Forwards Prediction Error can be used to choose lambda

# We need to define a matrix where each row gives the index of the observation
# time to start solving the ODE from, and the observation time to stop. 
whichtimes = cbind(1:31,11:41)

# then we call

FPE = forward.prediction.error(FhNtimes,FhNdata,coefs2,lik,proc,pars2,whichtimes)

# Usually we would do this over a selection of lambda values

lambdas = 10^seq(1:8)
FPEs = 0*lambdas
for(ilam in 1:length(lambdas)){
  t.Ores = Profile.LS(fhn.fun,FhNdata,FhNtimes,FhNpars,coefs0,FhNbasis,lambdas[ilam])
  FPEs[ilam] = forward.prediction.error(FhNtimes,FhNdata,t.Ores$coefs,lik,proc,t.Ores$pars,whichtimes)
}

# And look at these values

FPEs


# FPE makes use of a function IntegrateForward which interfaces a proc object to
# lsoda. We can use it directly to solve the ODE at the estimated parameters and compre
# this to data as follows (assuming we know x0):

x0 = c(-1,1)
names(x0) = FhNvarnames

sol1 = IntegrateForward(x0,FhNtimes,Ores2$pars,proc)
par(mar=c(5,5,1,1))
matplot(FhNtimes,FhNdata,pch=c('V','R'),cex=1.5,cex.lab=2.5,cex.axis=2.5)
matplot(sol1$times,sol1$states,type='l',lwd=3,add=TRUE)



# Some diagnostic plots help to work out if there are problems with the model fit to the data
# or with the smooth fit to the model. 


out1 = CollocInferPlots(Ores2$coefs,Ores2$pars,lik,proc,times=FhNtimes,data=FhNdata,
                            cex.lab=2.5,cex.axis=2.5,cex=1.5,lwd=3)


# Covariance of parameters can be obtained from 

covar = Profile.covariance(Ores2$pars,times=FhNtimes,data=FhNdata,coefs=Ores2$coefs,lik=lik,proc=proc)

covar

# and we can look at confidence intervals

CIs = cbind( Ores2$pars - 2*sqrt(diag(covar)), Ores2$pars + 2*sqrt(diag(covar)) )
rownames(CIs) = FhNparnames
CIs

## When only some state variables are observed, we can still conduct profiling. 
# to mimic this, first set the second column of data to be NA

data2 = FhNdata
data2[,2] = NA

# and remember that we also don't get coefficients, we'll set these to zero

coefs0.2 = coefs0
coefs0.2[,2] = 0

# The function FitMatchOpt allows us to pull some columns of the coefficient
# matrix into line with the differential equation (keeping the other columns fixed):

Fres3 = FitMatchOpt(coefs0.2,2,pars1,proc)

# And we can now proceed with profiling as before:

Ores4 = outeropt(FhNdata,FhNtimes,pars1,Fres3$coefs,lik,proc)
Ores4$pars

out2 = CollocInferPlots(Ores4$coefs,Ores4$pars,lik,proc,times=FhNtimes,data=data2)

covar4 = Profile.covariance(Ores4$pars,times=FhNtimes,data=data2,coefs=Ores4$coefs,lik=lik,proc=proc)
CI4 = cbind( Ores4$pars - 2*sqrt(diag(covar4)), Ores4$pars + 2*sqrt(diag(covar4)) )
rownames(CI4) = FhNparnames
CI4
