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
# The system has 16 parameters. Notable features include that only the sums 
# C1+C2 and B+S can be observed. Further, an unknown fraction of each is counted 
# at each time. This requires us to set up a model for the observation process 
# along with the ODE. 

# First we load up some data

data(ChemoData)

# The first two of these parameters give the fractions of Algae and Rotifers 
# that are counted. The remaining parameters are all positive and using 
# their logged values is helpful. 

logpars=c(ChemoPars[1:2],log(ChemoPars[3:16]))

# Parameters 'p1' and 'p2' represent relative palatability of the two algal
# clones, as such only one can be estimated and we fix p2 = 0. 

active = c(1:2,5,7:16)     


# We need some basis functions

rr = range(ChemoTime)
knots = seq(rr[1],rr[2],by=0.5)
bbasis = create.bspline.basis(rr,norder=4,breaks=knots)

# smooth obtained from generating some ODE solutions

y0 = log(c(2,0.1,0.4,0.2,0.1))  
names(y0) = ChemoVarnames

odetraj = lsoda(y0,ChemoTime,func=chemo.ode,logpars)

DEfd   = smooth.basis(ChemoTime,odetraj[,2:6],fdPar(bbasis,int2Lfd(2),1e-6))
coefs0 = DEfd$fd$coef


# We need a measurement model

ChemoMeas = function(t,x,p,more)
{
  return( cbind( log(exp(x[,'C1'])+exp(x[,'C2'])) + log(p['a1']), log(exp(x[,'B'])+exp(x[,'S']))+log(p['a2']) ) )
}


# From here we can employ the usual setup functions

out = LS.setup(pars = logpars,coefs = coefs0,fn=chemo.fun,
        basisvals=bbasis,lambda=c(1e3,2e2,2e2,1e3,1e3)/2,data=ChemoData,times=ChemoTime,
        posproc=TRUE,names=ChemoVarnames,likfn=ChemoMeas)


lik = out$lik
proc = out$proc


# Now, with parameters fixed, we'll estiamte coefficients. 

control.in = list()
control.in$trace = 2
control.in$maxit = 1000
control.in$reltol = 1e-6

res1 = inneropt(coefs=coefs0, pars=logpars, times=ChemoTime, data=log(ChemoData),
               lik=lik, proc=proc, in.meth='optim', control.in=control.in)

# We'll plot agreement with the data

out1 = CollocInferPlots(res1$coefs,logpars,lik,proc,ChemoTime,log(ChemoData))


# And conduct the outer optimization

res2 = outeropt(coefs=coefs0, pars=logpars, times=ChemoTime, data=log(ChemoData),
               lik=lik, proc=proc, in.meth='optim', control.in=control.in)
               
# along with diagnostic plots

out2 = CollocInferPlots(res2$coefs,res2$pars,lik,proc,ChemoTime,log(ChemoData))

###############################################################

# The following provides a manual setup that employs the builtin genin
# functions for measurement or process models that involve linear combinations
# of variables. 

mids = c(min(knots),(knots[1:(length(knots)-1)] + 0.25),max(knots))

bvals.obs = eval.basis(ChemoTime,bbasis)

bvals.proc = list(bvals  = eval.basis(mids,bbasis),
                  dbvals = eval.basis(mids,bbasis,1));


# We'll choose a fairly large value of lambda. 

lambda = c(1e3,2e2,2e2,1e3,1e3)


# We can now set up the proc object. We will want to take a log transformation
# of the state here for numerical stability. In general it is better to do 
# finite differencing AFTER the log transformation rather than before it. 

proc       = make.SSEproc()              # Sum of squared errors
proc$bvals = bvals.proc                  # Basis values

proc$more          = make.findif.ode()   # Finite differencing
proc$more$qpts     = mids                # Quadrature points
proc$more$weights  = rep(1,5)*lambda     # Quadrature weights (including lambda)
proc$more$names    = ChemoVarnames       # Variable names
proc$more$parnames = ChemoParnames       # Parameter names

proc$more$more = list(fn=make.logtrans()$fn,eps=1e-8) # Log transform

proc$more$more$more = list(fn=chemo.fun) # ODE function

# For the lik object we need to both represent the linear combination transform
# and we need to model the observation process. 

# First to represent the observation process, we can use the genlin
# functions. These produce a linear combination of the the states 
# (they can be used in proc objects for linear systems, too).  

temp.lik      = make.SSElik()
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
#temp.lik$more$weights = matrix(c(100,1),length(ChemoTime),2)

temp.lik$more$weights = c(100,1)

# Finally, we tell CollocInfer that the trajectories are represented on
# the log scale and must be exponentiated before comparing them to the data. 

lik       = make.logstate.lik()
lik$more  = temp.lik
lik$bvals = bvals.obs

# Now lets try running this

# Because we don't have direct observations of any state, we'll use a starting 
# smooth obtained from generating some ODE solutions

y0 = log(c(2,0.1,0.4,0.2,0.1))  
names(y0) = ChemoVarnames

odetraj = lsoda(y0,ChemoTime,func=chemo.ode,logpars)

DEfd   = smooth.basis(ChemoTime,odetraj[,2:6],fdPar(bbasis,int2Lfd(2),1e-6))
coefs0 = DEfd$fd$coef

# Now, with parameters fixed, we'll estiamte coefficients. 

control.in = list()
control.in$trace = 2
control.in$maxit = 1000
control.in$reltol = 1e-6

res = inneropt(coefs=coefs0, pars=logpars, times=ChemoTime, data=ChemoData,
               lik=lik, proc=proc, in.meth='optim', control.in=control.in)

# We'll for the trajectory and also the appropriate sum of exponentiated
# states to compare to the data. 

coefs1  = matrix(res$coefs,dim(coefs0))
traj    = lik$bvals %*% coefs1
obstraj = lik$more$more$fn(ChemoTime,exp(traj),logpars,lik$more$more$more)

# Plot these against the data

X11()
par(mfrow=c(2,1))
plot(obstraj[,1],type='l',ylab='Chlamy',xlab='',cex.lab=1.5,cex.axis=1.5)
points(ChemoData[,1])
plot(obstraj[,2],type='l',ylab='Brachionus',xlab='days',cex.lab=1.5,cex.axis=1.5)
points(ChemoData[,2])


# Now we can continue with the outer optimization

res2 = outeropt(pars=logpars,times=ChemoTime,data=ChemoData,coef=coefs1,
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