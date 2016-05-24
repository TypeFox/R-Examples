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


# First we generate some data 

data(ChemoData)

# The first two of these parameters give the fractions of Algae and Rotifers 
# that are counted. The remaining parameters are all positive and using 
# their logged values is helpful. 

logpars=c(ChemoPars[1:2],log(ChemoPars[3:16]))

# Parameters 'p1' and 'p2' represent relative palatability of the two algal
# clones, as such only one can be estimated and we fix p2 = 0. 

active = c(1:2,5,7:16)     

# We'll choose a fairly large value of lambda. 

lambda = rep(10000,5)  

# We need some basis functions

rr = range(ChemoTime)
knots = seq(rr[1],rr[2],by=0.5)
bbasis = create.bspline.basis(rr,norder=4,breaks=knots)

# We now need to define two functions. First, chemo.fun defines a map from the
# trajectory to its derivative

chemo.fun = function (times, y, p, more = NULL) 
{
    r = y
    p = exp(p)
    Q = p["p1"] * y[, "C1"] + p["p2"] * y[, "C2"]
    Qs = Q * exp(100 * (Q - p["Qstar"]))/(1 + exp(100 * (Q - 
        p["Qstar"]))) + p["Qstar"]/(1 + exp(100 * (Q - p["Qstar"])))
    r[, "N"] = p["delta"] * (p["NI"] - y[, "N"]) - p["rho"] * 
        y[, "C1"] * y[, "N"]/(p["KC1"] + y[, "N"]) - p["rho"] * 
        y[, "C2"] * y[, "N"]/(p["KC2"] + y[, "N"])
    r[, "C1"] = y[, "C1"] * (p["XC"] * p["rho"] * y[, "N"]/(p["KC1"] + 
        y[, "N"]) - p["p1"] * p["G"] * (y[, "B"] + y[, "S"])/(p["KB"] + 
        Qs) - p["delta"])
    r[, "C2"] = y[, "C2"] * (p["XC"] * p["rho"] * y[, "N"]/(p["KC2"] + 
        y[, "N"]) - p["p2"] * p["G"] * (y[, "B"] + y[, "S"])/(p["KB"] + 
        Qs) - p["delta"])
    r[, "B"] = y[, "B"] * (p["XB"] * p["G"] * Q/(p["KB"] + Qs) - 
        (p["delta"] + p["m"] + p["lambda"]))
    r[, "S"] = p["lambda"] * y[, "B"] - (p["delta"] + p["m"]) * 
        y[, "S"]
    return(r)
}

# Second, we need a map from the trajectory to the observations. We'll use 
# log observations in this case since the observations are positive and this
# helps to stabilize scales. 

chemo.obs = function(t,y,p,more)
{
  Algae = p[1]*(y[,'C1']+y[,'C2'])
  Rotifers = p[2]*(y[,'B']+y[,'S'])
  
  return( log( cbind(Algae,Rotifers) ) )

}


### From here we will set-up the objects needed to run profiling. Below we give
# LS.setup the parameters, maps from trajectories to derivatives and to
# observations, along with data, times, basis functions, variable names and
# trade-off parameter lambda. 
#
# Additionally, posproc and poslik indicate that we are modelling the log of the
# trajectory, so that we should exponentiate before mapping to derivatives and
# before mapping to observations. 

out = LS.setup(pars = logpars, fn = chemo.fun, basisvals=bbasis,lambda=lambda,
               data = log(ChemoData),times=ChemoTime,posproc=TRUE,poslik=TRUE,
               names = ChemoVarnames,likfn=chemo.obs)

lik = out$lik
proc = out$proc

# Because we don't have direct observations of any state, we'll use a starting
# smooth obtained from generating some ODE solutions. 

y0 = log(c(2,0.1,0.4,0.2,0.1))
names(y0) = ChemoVarnames

odetraj = IntegrateForward(y0,ChemoTime,logpars,proc)

DEfd   = smooth.basis(odetraj$times,odetraj$states,fdPar(bbasis,int2Lfd(2),1e-6))
coefs0 = DEfd$fd$coef

obsvals = lik$more$fn(ChemoTime,odetraj$states,ChemoPars,lik$more$more)

ChemoData  = exp(obsvals + 0.5*matrix(rnorm(2*length(ChemoTime)),length(ChemoTime),2))



## We can now call the inner optimization -- this just tries to fined the 
# coefficients that both fit the data and the differential equation well. 

res = inneropt(coefs=coefs0, pars=logpars, times=ChemoTime, data=log(ChemoData),
               lik=lik, proc=proc, in.meth='optim',control.in=list(trace=2))

# We'll form the trajectory and also the appropriate sum of exponentiated
# states to compare to the data. 

traj    = lik$bvals %*% res$coefs
obstraj = lik$more$fn(ChemoTime,traj,logpars,lik$more$more)

# Plot these against the data

X11()
par(mfrow=c(2,1))
plot(obstraj[,1],type='l',ylab='Chlamy',xlab='',cex.lab=1.5,cex.axis=1.5)
points(log(ChemoData[,1]))
plot(obstraj[,2],type='l',ylab='Brachionus',xlab='days',cex.lab=1.5,cex.axis=1.5)
points(log(ChemoData[,2]))

# The function outeropt now allows parameter estimation to take place through 
# profiling procedure. 

res2 = outeropt(pars=logpars,times=ChemoTime,data=log(ChemoData),coef=res$coefs,
    lik=lik,proc=proc,active=active,in.meth='optim',out.meth='nlminb')

# We'll extract the resulting parameters and coefficients. 

npars = res2$pars
C = res2$coefs

# And obtain an estimated trajectory and the exponentiated sum to comprare
# to the data. 

traj = lik$bvals%*%C
ptraj = lik$more$fn(ChemoTime,traj,npars,lik$more$more)

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
points(ChemoTime,log(ChemoData[,1]))
plot(ChemoTime,ptraj[,2],type='l',ylab='Brachionus',xlab='days',cex.lab=1.5)
points(ChemoTime,log(ChemoData[,2]))


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
  plot(proc$more$qpts,dtraj2[,i],type='l',xlab='',ylab=ChemoVarnames[i],
    cex.lab=1.5,ylim=c(-0.5,0.5))
  lines(proc$more$qpts,ftraj2[,i],col=2,lty=2)
  abline(h=0)
}
legend('topleft',legend=c('Smooth','Model'),lty=1:2,col=1:2,cex=1.2)

# Solving the differential equation from the estiamted initial condition
# of the trajectory allows us to compare the qualitative behavior of 
# our estimate to that of the differential equation. 

y0 = traj[1,]
names(y0) = ChemoVarnames
odetraj = IntegrateForward(y0,ChemoTime,npars,proc)


X11()
par(mfrow=c(2,1))
matplot(ChemoTime,traj,col=1,type='l',lwd=3,cex.lab=1.5,cex.axis=1.5,
    ylab='',cex.main=1.5,main='Reconstructed Trajectories')
legend(x='topright',legend=ChemoVarnames,lwd=3,col=1,lty=1:5)
matplot(ChemoTime,odetraj$states,col=1,type='l',lwd=3,cex.axis=1.5,cex.lab=1.5,
    ylab='',cex.main=1.5,main='ODE Solution')

# We can also compare the pattern of observations predicted by the differential
# equation and that estimated by our methods. 

otraj = lik$more$fn(ChemoTime,odetraj$states,npars,lik$more$more)

X11()
par(mfrow=c(2,1))
matplot(ChemoTime,ptraj,type='l',lwd=2,xlab='days',cex.lab=1.5,ylab='',
  cex.axis=1.5,cex.main=1.5,main='Predicted Observations -- Smooth')
matplot(ChemoTime,log(ChemoData),add=TRUE,pch = c(1,2))
legend('topright',legend=c('Algae','Rotifers'),pch=1:2,col=1:2)
matplot(ChemoTime,otraj,type='l',lwd=2,xlab='days',cex.lab=1.5,ylab='',
  cex.axis=1.5,cex.main=1.5,main='Predicted Observations -- ODE')
legend('topright',legend=c('Algae','Rotifers'),lty=1:2,col=1:2,lwd=2)