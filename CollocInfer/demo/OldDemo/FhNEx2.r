# This demonstration file is intended to accompany the CollocInfer manual. 

library(CollocInfer)

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

fhn.ode <- function(times,y,p)
{
 r = y;
 dimnames(r) = dimnames(y);
 r['V'] = p['c']*(y['V'] - y['V']^3/3 + y['R'])
 r['R'] = -(y['V'] -p['a'] + p['b']*y['R'])/p['c']

 return(list(r))
}
 
# We need the times at which we will observe the system
 
FhNtimes = seq(0,20,0.05)

# And now we can create solutions to the equations

out = lsoda(x0,times=FhNtimes,fhn.ode,FhNpars)

matplot(out[,1],out[,2:3],type='l',lwd=2)
legend('bottomleft',c('V','R'),lwd=2,col=1:2,lty=1:2)

plot(out[,2],out[,3],type='l',lwd=2,xlab='V',ylab='R')


FhNdata = out[,2:3] + 0.06*matrix(rnorm(802),length(FhNtimes),2)

rm(list=ls())
data(FhNdata)

# Define a basis

range = c(0,20)
knots = seq(0,20,0.5)
norder = 4

FhNbasis = create.bspline.basis(range=range,norder=norder,breaks=knots)

# Smooth out the data

DEfd0 = smooth.basis(FhNtimes,FhNdata,fdPar(FhNbasis,int2Lfd(2),1))$fd

# Produce a plot of how well this agrees with the smooth
par(mfrow=c(2,1))
plotfit.fd(FhNdata, FhNtimes, DEfd0)


coefs0 = DEfd0$coef
colnames(coefs0) = FhNvarnames

# Define right hand side function

fhn.fun <- function(times,y,p,more)
{
 r = y;
 r[,'V'] = p['c']*(y[,'V'] - y[,'V']^3/3 + y[,'R'])
 r[,'R'] = -(y[,'V'] -p['a'] + p['b']*y[,'R'])/p['c']
 return(r)
}
     
# Now choose a lambda and obtain profiling values     


lambda = 1000

profile.obj = LS.setup(pars=FhNpars,fn=fhn.fun,lambda=lambda,times=FhNtimes,
                       coefs=coefs0,basisvals=FhNbasis)
                       
proc = profile.obj$proc
lik = profile.obj$lik
                       
# Start with gradient matching

Pres0  = ParsMatchOpt(FhNpars,coefs0,proc)
pars1 = Pres0$pars

# Now move on to inner optimization

Ires1 = inneropt(FhNdata,times=FhNtimes,pars1,coefs0,lik,proc)
coefs1 = Ires1$coefs

# Alternatively


Ires1.2 = Smooth.LS(fhn.fun,FhNdata,FhNtimes,FhNpars,coefs0,FhNbasis,lambda)

# And outer optimization

Ores2 = outeropt(FhNdata,FhNtimes,pars1,coefs1,lik,proc)

Ores2$pars


# Alternatively

Ores2.2 = Profile.LS(fhn.fun,FhNdata,FhNtimes,FhNpars,coefs0,FhNbasis,lambda)


# Forwards Prediction Error

whichtimes = cbind(1:31,11:41)

FPE = forward.prediction.error(FhNtimes,FhNdata,Ores2$coefs,lik,proc,Ores2$pars,whichtimes)

lambdas = 10^seq(1:8)
FPEs = 0*lambdas
for(ilam in 1:length(lambdas)){
  t.Ores = Profile.LS(fhn.fun,FhNdata,FhNtimes,FhNpars,coefs0,FhNbasis,lambdas[ilam])
  FPEs[ilam] = forward.prediction.error(FhNtimes,FhNdata,t.Ores$coefs,
                                   lik,proc,t.Ores$pars,whichtimes)
}


# Also IntegrateForward

x0 = c(-1,1)
names(x0) = FhNvarnames

sol1 = IntegrateForward(x0,FhNtimes,Ores2$pars,proc)
matplot(FhNtimes,FhNdata,cex.lab=1.5,cex.axis=1.5)
matplot(sol1$times,sol1$states,type='l',add=TRUE)

# Some diagnostic plots

out1 = CollocInferPlots(Ores2$coefs,Ores2$pars,lik,proc,times=FhNtimes,data=FhNdata)


# Covariance of parameters

covar = Profile.covariance(Ores2$pars,times=FhNtimes,data=FhNdata,coefs=Ores2$coefs,
                             lik=lik,proc=proc)

# And obtain confidence intervals

CIs = cbind( Ores2$pars - 2*sqrt(diag(covar)), Ores2$pars + 2*sqrt(diag(covar)) )
rownames(CIs) = FhNparnames
CIs

## When only some state variables are observed:

data2 = FhNdata
data2[,2] = NA

coefs0.2 = coefs0
coefs0.2[,2] = 0

Fres3 = FitMatchOpt(coefs0.2,2,pars1,proc)

Ores4 = outeropt(FhNdata,FhNtimes,pars1,Fres3$coefs,lik,proc)
Ores4$pars

out2 = CollocInferPlots(Ores4$coefs,Ores4$pars,lik,proc,times=FhNtimes,data=data2)

covar4 = Profile.covariance(Ores4$pars,times=FhNtimes,data=FhNdata,
                                    coefs=Ores4$coefs,lik=lik,proc=proc)
                                    
CI4 = cbind( Ores4$pars - 2*sqrt(diag(covar4)), Ores4$pars + 2*sqrt(diag(covar4)) )