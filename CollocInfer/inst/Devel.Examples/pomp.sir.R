require(pomp)
library(CollocInfer)


## euler.sim pomp object from the pomp library

data(euler.sir)

## Parameters

pars = coef(euler.sir)
parnames = names(pars)


## Observations and observation times; we'll start by assuming that the 
## observations are given by a scaled version of the infectives. 

Y = t(data.array(euler.sir))
times = time(euler.sir)

Data = matrix(NA,nrow(Y),5)   ## Easiest way to specify that we only observe
Data[,2] = log(Y) + 1.217     ## the infectives. 


## Solve the ODE to get a trajectory and state names

x <- trajectory(euler.sir)
varnames =  dimnames(x)[[1]]  ## I really wish I could extract this sort of thing
x = t(x[,1,])                 ## directly from the pomp object

x[,5] = 1

## Now set up a basis representation for the ODE and start by smoothing the 
## solved trajectory

ebasis = create.bspline.basis(c(0,4),norder=4,breaks=seq(0,4,0.1))

trajfd = smooth.basis(times,log(x[-1,]),fdPar(ebasis,int2Lfd(2),0.01))

coefs = trajfd$fd$coefs
colnames(coefs) = varnames

## Now we'll set up the CollocInfer objects

objs = LS.setup(pars=pars, coefs=coefs, fn=euler.sir,
                basisvals=ebasis, lambda = 1e1, 
                data = Data, times=times, posproc=TRUE)

## And try smoothing

res1 = inneropt(Data,times,pars,coefs,objs$lik,objs$proc,control.in=list(trace=6))


## And the outer optimization

res2 = outeropt(Data,times,pars,res1$coefs,objs$lik,objs$proc,active=1:3)

##########################
## Try without logging
#############################

Data = exp(Data)

trajfd = smooth.basis(times,x[-1,],fdPar(ebasis,int2Lfd(2),0.00001))

coefs = trajfd$fd$coefs
colnames(coefs) = varnames

coefs[,2] = smooth.basis(times,Data[,2],fdPar(ebasis,int2Lfd(2),0.00001))$fd$coef


## Now we'll set up the CollocInfer objects

objs = LS.setup(pars=pars, coefs=coefs, fn=euler.sir,
                basisvals=ebasis, lambda = 1e-3, 
                data = Data, times=times, posproc=FALSE)

## And try smoothing

res1 = inneropt(Data,times,pars,coefs,objs$lik,objs$proc,control.in=list(trace=6))

## And the outer optimization

res2 = outeropt(Data,times,pars,res1$coefs,objs$lik,objs$proc,active=1:3)



##### Now something more curious on the observation process; we'll use the difference
##### in the cases state variable.

proc = objs$proc
lik = objs$lik

lik$bvals = eval.basis(times,ebasis) - eval.basis(c(0,times[-208]),ebasis)

Data = matrix(0,208,5)
Data[,4] = Y

## And try smoothing

res1 = inneropt(Data,times,pars,coefs,lik,proc,control.in=list(trace=6))

## And the outer optimization

res2 = outeropt(Data,times,pars,res1$coefs,lik,proc,active=1:3)
