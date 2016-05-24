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



