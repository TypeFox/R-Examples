## Example Diagnostics -- Learning the FitzHugh-Nagumo Equations
#
# Demonstration estimation of function functions. 

library('CollocInfer')

# First Create Some data from the FitzHugh-Nagumo equations. We will first estimate
# a linear differential equation to mimic this, and then we will try to work 
# out how to fix it up

t = seq(0,20,0.05)

pars = c(0.2,0.2,3)
names(pars) = c('a','b','c')

x0 = c(-1,1)
names(x0)= c('V','R')
y = lsoda(x0,t,make.fhn()$fn.ode,pars)
y = y[,2:3]

data = y + 0.2*matrix(rnorm(802),401,2)

# Now define a basis object

knots = seq(0,20,0.2)
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	norder=norder,breaks=knots)


# Initial values for coefficients will be obtained by smoothing

DEfd = smooth.basis(t,data,fdPar(bbasis,1,0.5)) 
coefs = DEfd$fd$coefs

names(coefs) = c('V','R')

# Starting parameter estimates correspond to circular motion

spars = c(0,1,-1,0)

# Now set up some profiling; make.genlin() produces a linear ODE using parameters
# in row-wise order in the equaiton Dx = Ax

out = LS.setup(coefs=coefs,pars=spars,times=t,fn=make.genlin(),basisvals=bbasis,
                  lambda=c(100,100),names=c('V','R'))

lik = out$lik
proc = out$proc

# With this we can run model-based smoothing

Ires = inneropt(data,times=t,spars,coefs,lik,proc)

# and profiling

Ores = outeropt(data=data,times=t,pars=spars,coefs=Ires$coefs,lik=lik,proc=proc,
  in.meth="nlminb",out.meth="nlminb")

# There is some clear lack of fit here

out2 = CollocInferPlots(Ores$coefs,Ores$pars,lik,proc,times=t,data=data)

# And we will attempt to relate the mis-match in Dx and f(x) to the value of x
# to see how we might be able to fix them up. 

par(mfrow=c(2,2))
plot( out2$traj[,1], out2$dtraj[,1]-out2$ftraj[,1],type='l')
plot( out2$traj[,2], out2$dtraj[,1]-out2$ftraj[,1],type='l')
plot( out2$traj[,1], out2$dtraj[,2]-out2$ftraj[,2],type='l')
plot( out2$traj[,2], out2$dtraj[,2]-out2$ftraj[,2],type='l')

## Now set up empirical forcing functions to get a more precise handle on this
# by examining the equations  Dx = f(x) + Psi(t) P where f(x) is the equation we
# estimated above, Psi is a new basis system and we will estimate the coefficients
# P as parameters via the same method.               

# The additional functions g(t) = Psi(t) P are known as "empirical forcing
# functions", hence the "f" in front of various objects employed below.

# First we need a set of basis functions at which to evaluate the forcing
# functions and their evaluation at the quadrature points

fbasis = create.bspline.basis(range=range,nbasis=23,norder=4)
fbvals = eval.basis(proc$more$qpts,fbasis)

# Now we can call the steup functions, in this case we need to give a lot
# to the more object. In particular, we start off with the first more
# element of the proc object we just used

procmore = proc$more

# and we have to add to this the estimated parameters

procmore$p = Ores$pars

# and which elements to add the forcing functions to

procmore$which = 1:2

# and the evalution of fbasis at the quadrature points

procmore$psi = fbvals

# Now we call the setup function

out3 = LS.setup(coefs=Ores$coefs,pars=rep(0,fbasis$nbasis),times=t,
                  fn=make.diagnostics(),basisvals=bbasis,
                  lambda=c(100,100),names=c('V','R'),more=procmore)
                  
lik2 = out3$lik
proc2 = out3$proc                  


# With this we can run model-based smoothing

Ires2 = inneropt(data,times=t,rep(0,2*fbasis$nbasis),Ores$coefs,lik2,proc2)

# and profiling

Ores2 = outeropt(data=data,times=t,pars=rep(0,2*fbasis$nbasis),coefs=Ires2$coefs,
            lik=lik2,proc=proc2,in.meth="nlminb",out.meth="nlminb")

# Usual CollocInferPlots

out4 = CollocInferPlots(Ores2$coefs,Ores2$pars,lik2,proc2,times=t,data=data)

# Now we can look at the values of the time-varying forcing functions. These 
# values are given by combining the forcing function basis and the parameters

fvals = eval.basis(out4$timevec,fbasis)%*%matrix(Ores2$pars,fbasis$nbasis,2)

# which we can plot

matplot(out4$timevec,fvals,type='l')

# And relate to the estimated trajectories where the mising cubic term can
# be discerned.

par(mfrow=c(2,1))
matplot(out4$traj[,1],fvals,type='l')
matplot(out4$traj[,2],fvals,type='l')




## Manual set-up. 

# This is the equivalent analysis to what we have done above, but
# we create the lik and proc objects manually -- we believe that the structures
# we work with should not be hidden from the user. 

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



