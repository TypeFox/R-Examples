# This demonstration file is intended to accompany the CollocInfer manual. 

library(CollocInfer)

# To ensure reproducibility
set.seed((2004*2007)/2014)

# This file carrys out the analysis for the Rosenzweig-MacArthur model that is
# given in the manual. We assume that the reader is familiar with the FitzHugh-Nagumo
# example that is first developed in the manual, so we will focus on the elaborations
# that this model requires. 

# First, we write down the three-species Rosenzweig-Macarthur model in a form 
# suitable for CollocInfer

RosMac2 = function(t,x,p,more){

  p = exp(p)
  dx = x

  dx[,'C1'] = p['rho1']*x[,'C1']*(1- x[,'C1']/p['kappaC1']- x[,'C2']/p['kappaC2']) - p['pi']*p['gamma']*x[,'C1']*x[,'B']/(p['kappaB']+p['pi']*x[,'C1']+x[,'C2'])
  dx[,'C2'] = p['rho2']*x[,'C2']*(1- x[,'C1']/p['kappaC1']- x[,'C2']/p['kappaC2']) - p['gamma']*x[,'C2']*x[,'B']/(p['kappaB']+p['pi']*x[,'C1']+x[,'C2'])
  dx[,'B'] = p['chi']*p['gamma']*(p['pi']*x[,'C1']+x[,'C2'])*x[,'B']/(p['kappaB']+p['pi']*x[,'C1']+x[,'C2']) - p['delta']*x[,'B']

  return(dx)
}

# We will define some interesting parameters

RMpars = c(0.2,0.025,0.125,2.2e4,1e5,5e6,1,1e9,0.3)
RMParnames = c('pi','rho1','rho2','kappaC1','kappaC2','gamma','chi','kappaB','delta') 

# Which we represent on the log scale

logpars = log(RMpars)
names(logpars) = RMParnames

# And we also need some initial conditions (with named entries)

RMVarnames = c('C1','C2','B')

x0 = c(50,50,2)
names(x0) = RMVarnames

# The following version of RosMac2 is suitable for use with lsoda

RosMac2ODE = function(t,z,p){
  p = exp(p)
  x = exp(z)
  dx = x

  dx['C1'] = p['rho1']*x['C1']*(1- x['C1']/p['kappaC1']-x['C2']/p['kappaC2']) - p['pi']*p['gamma']*x['C1']*x['B']/(p['kappaB']+p['pi']*x['C1']+x['C2'])
  dx['C2'] = p['rho2']*x['C2']*(1- x['C2']/p['kappaC2']- x['C1']/p['kappaC1']) - p['gamma']*x['C2']*x['B']/(p['kappaB']+p['pi']*x['C1']+x['C2'])
  dx['B'] = p['chi']*p['gamma']*(p['pi']*x['C1']+x['C2'])*x['B']/(p['kappaB']+p['pi']*x['C1']+x['C2']) - p['delta']*x['B']

  return(list(dx/x))
}

# With this we can solve the ODE at 200 successive days
time = 0:200
res0 = lsoda(log(x0),time,RosMac2ODE,p = logpars)

# and plot the solutions
matplot(exp(res0[,2:4]),type='l')

# We'll obtain data by adding noise
data = res0[,2:4] + 0.2*matrix(rnorm(603),201,3)

# and name the columns
colnames(data) = RMVarnames

# Giving us the following plot
matplot(data,cex.lab=2.5,cex.axis=2.5,cex=1.5,xlab='days',pch = c('1','2','B'))
matplot(res0[,2:4],type='l',lwd=3,add=TRUE)


# Now we need to set up the CollocInfer machinery

# First we'll define a basis with knots each time point

rr = range(time)
knots = seq(rr[1],rr[2],by=1)

bbasis = create.bspline.basis(rr,norder=4,breaks=knots)

# And obtain an initial set of parameters and coefficients from smoothing

coef0 = smooth.basis(time,data,fdPar(bbasis,int2Lfd(2),10))$fd$coef
colnames(coef0) = RMVarnames

# We will now create the profiling objects, but use the log transformation by
# setting posproc=TRUE 

out = LS.setup(pars=logpars,coefs=coef0,basisvals=bbasis,fn=RosMac2,lambda=1e5,
          times=time,posproc=TRUE)
          
          
lik = out$lik
proc = out$proc

# We'll do gradient matching to get parameter estimates corresponding to this
# smooth

res1 = ParsMatchOpt(logpars,coef0,proc)

# And now profiling
res3 = outeropt(data,time,res1$pars,coef0,lik,proc)

# Let's have a look at the parameters that we got
exp(res3$pars)

# and the agreement with data and model
out3 = CollocInferPlots(res3$coefs,res3$pars,lik,proc,times=time,data=data,
                  cex.lab=2.5,cex.axis=2.5,cex=1.5,lwd=3)


## Now we will compliate the model. In reality, we only observe the sum of C1 
# and C2. 

data2 = cbind( log( exp(data[,'C1'])+exp(data[,'C2'])), data[,'B'])
matplot(data2,cex.lab=2.5,cex.axis=2.5,pch=c('C','B'),cex=1.5,xlab='days')

# To deal with this, we need to define a transformation function from our 
# state variables to the expected observations. In this case (since the states
# are on the log scale) we exponentiate to get back to the original scale, add
# C1 and C2 and then take the log again. 
RMobsfn = function(t,x,p,more)
{
  x = exp(x)
  y = cbind( x[,'C1']+x[,'C2'],x[,'B'])
  return(log(y))
}

# We can now create new profiling objects that incorporate this transformation
# function by specfying likfn

out = LS.setup(pars=logpars,coefs=coef0,basisvals=bbasis,fn=RosMac2,lambda=1e5,
          times=time,posproc=TRUE,likfn=RMobsfn)
          
lik2 = out$lik
proc2 = out$proc

# To see how we perform in this situation, we'll start by setting the two 
# columns of the data smooth to zero

coef02 = coef0
coef02[,1:2] = 0

# We'll now pull these columns into line with the rotifer column (which we can
# still smooth. 

Fres3 = FitMatchOpt(coef02,1:2,res1$pars,proc2)

# And we'll run profiling and have a look at what we get. 
res32 = outeropt(data2,time,res1$pars,Fres3$coefs,lik2,proc2)
exp(res32$pars)
out32 = CollocInferPlots(res32$coefs,res32$pars,lik2,proc2,times=time,data=data2,
                  datanames=c('C','B'),cex.lab=2.5,cex.axis=2.5,cex=1.5,lwd=3)

# The section at the end of this demo goes through setting up lik2 and proc2
# manually rather than through LS.setup. 


### In this framework is is not unreasonable to expect that we have
# repeated experiments. When these are very well structured and all have common
# observation times, this can be easily accommodated in CollocInfer (it can 
# accommodate less regular replicated experiments, but requires work to set things
# up manuall. 

# We'll create a second experiment starting from new initial conditions

x03 = c(15,25,4)
names(x03) = RMVarnames

res03 = lsoda(log(x03),time,RosMac2ODE,p = logpars)

data03 =  res03[,2:4] + 0.2*matrix(rnorm(603),201,3)

# and set up a three dimensional array in which the experiment number is the
# second dimension and the third dimension is the variable being measured (this
# is to agree with conventions in the fda package)

alldat = array(0,c(201,2,3))
alldat[,1,] = data
alldat[,2,] = data03

# Nowe we'll smooth the second experiment

coef3 = smooth.basis(time,data03,fdPar(bbasis,int2Lfd(2),10))$fd$coef

# and create a three-dimensional array with all the coefficients together

coefs = array(0,c(dim(coef3)[1],2,3))
coefs[,1,] = coef0
coefs[,2,] = coef3

# These three dimensional arrays can be given to LSsetup which understands
# thi structure and knows what to do with it. 

out = LS.setup(pars=logpars,coefs=coefs,basisvals=bbasis,fn=RosMac2,lambda=1e5,
          times=time,data=alldat,posproc=TRUE,names=RMVarnames)
          
          
lik3 = out$lik
proc3 = out$proc

# We can now call the inner optimisation to use a model-based smooth
res13 = inneropt(data=out$data,times=out$times,pars=res1$pars,coefs=out$coefs,lik=lik3,proc=proc3)

# And use profiling to estimate parameters
res33 = outeropt(data=out$data,times=out$times,res1$pars,res13$coefs,lik3,proc3)

# These parameters should hopefully be closer to the truth than with only one experimental run
exp(res33$pars)

# And we can also examine the fit to the data and model. In this case, the times vector
# wraps around creating a few unpleasant graphical effects. 

out3 = CollocInferPlots(res33$coefs,res33$pars,lik3,proc3,times=out$times,data=out$data,datanames=c('B','C'),
                  cex.lab=2.5,cex.axis=2.5,cex=1.5,lwd=3)



###### Manual set-up

# Here we create lik2 and proc2 (corresponding to the indirectly observed 
# single-run experiment above in order to demonstrate their structure. 

# First we need to specify the matrices of basis values that we will use

# at observation time points
bvals.obs = eval.basis(time,bbasis)

# and quadrature times, these are midpoints between knots plus the end points
# The quadrature weights are all equal, but we have multiplied them by the lambda
# that we are using, in this case 1e5. 

qpts = c(knots[1],knots[1:(length(knots)-1)]+diff(knots)/2,knots[length(knots)])
qwts = 1e5*matrix(1,length(qpts),3)/length(qpts)

# basis values for proc is a list

bvals.quad = list(bvals  = eval.basis(qpts,bbasis), 
                  dbvals = eval.basis(qpts,bbasis,1))


# Now we create the lik object

# make.SSElik() sets up the squared error criteria
lik.m = make.SSElik()

# Attach the values of the basis expansion at the observation time points
lik.m$bvals = bvals.obs

# We need to specify the transformation of the state variables that is to be
# compared with the data. In this case, we will use finite differencing to 
# to compute the derivatives that we need; we can achieve this by first 
# employing make.findif.ode()
lik.m$more = make.findif.ode()

# and then telling findif.ode that the function it is finite differencing is
# RMobsfn
lik.m$more$more = list(fn = RMobsfn,eps=1e-6,more=NULL) 

# We also need to give lik.m a set of weights. This has to occur in 
# the more element of lik.m because it is used inside lik.m$fn. 
lik.m$more$weights = array(1,dim(data))


# We can also create the proc object manually. First we call make.SSEproc()
# in order to set up the squared error functions
proc.m = make.SSEproc()

# We also specify the basis values and their derivatives at the quadrature points
proc.m$bvals = bvals.quad

# Now we need to tell SSEproc about the right hand side of the ODE and its
# derivatives. Here we will use finite differencing again, as in lik.m:
proc.m$more = make.findif.ode()

# In fact, we are going to finite difference the right hand side of the 
# the ODE for the log transformed data. Here we specify the log transform
proc.m$more$more = list(fn = make.logtrans()$fn,eps=1e-6)

# and then give it the (non log-transform) Rosenzweig-MacArthur equations
proc.m$more$more$more$fn = RosMac2

# proc.m$more also needs to include some elements for internal processing. In 
# particular the following specify the quadrature points and weights
proc.m$more$weights = qwts
proc.m$more$qpts = qpts

# We will also tell it about the parameter and variable names. 
proc.m$more$parnames = RMParnames
proc.m$more$names = RMVarnames



# And let's check that this all works

Fres.m3 = FitMatchOpt(coef02,1:2,res1$pars,proc.m)
res.m32 = outeropt(data2,time,res1$pars,Fres3$coefs,lik2,proc.m)
exp(res.m32$pars)
out.m32 = CollocInferPlots(res.m32$coefs,res.m32$pars,lik.m,proc.m,times=time,data=data2)

# If you compare this to out32, you should have the same estimates. 
