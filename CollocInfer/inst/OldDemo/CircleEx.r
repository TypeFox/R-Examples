library('fda')
library('odesolve')
library('maxLik')
library('MASS')
library('Matrix')
library('SparseM')
source('../R/ProfileR.R')
source('../R/sse.shortcut.R')
source('Circle.R')
source('../R/findif.ode.R')
source('../R/SSElik.R')
source('../R/SSEproc.R')
source('../R/makeid.R')
source('../R/inneropt.R')
source('../R/cproc.shortcut.R')
source('../R/Cproc.R')
source('../R/genlin.R')
source('../R/cvar.R')
source('../R/Multinorm.R')

# First Create Some Data

t = seq(0,20,0.05)
pars = c(1,1)
y = array(0,c(length(t),2))
y[,1] =10*sin(t)
y[,2] =10*cos(t)
data = y + 0.05*matrix(rnorm(length(t)*2),length(t),2)

# Now a basis object

knots = seq(0,20,0.2)
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	norder=norder,breaks=knots)


# Initial values for coefficients will be obtained by smoothing

fd.data = array(data,c(dim(data)[1],1,dim(data)[2]))

DEfd = data2fd( fd.data,t,bbasis,fdnames=list(NULL,NULL,NULL) )

coefs = matrix(DEfd$coefs,dim(DEfd$coefs)[1],dim(DEfd$coefs)[3])
colnames(coefs) = DEfd$fdnames[[3]]



###Now lets try some optimization

spars =c(1,1.1)        # Perturbed parameters
lambda =10000
control=list()                # Control parameters 
control$trace = 0
control$maxit = 1000
control$maxtry = 10
control$reltol = 1e-6
control$meth = "BFGS"

control.in = control
control.in$reltol = 1e-12
control.out = control
control.out$trace = 2

control.in$print.level = 0
control.in$iterlim = 1000


#### set up the lik and proc objects ####
profile.obj = sse.setup(pars=pars,coefs,make.Circle(),basisvals=bbasis,lambda,fd.obj=NULL,more=NULL,
    weights=NULL,times=t,quadrature=NULL,eps=1e-6)
 lik = profile.obj$lik
proc = profile.obj$proc
coefs = profile.obj$coefs



res0 = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
        control=control.out,times=t,data=data,lik=lik,proc=proc,pars=spars)

ncoefs = matrix(res0$par,102,2)

### Multinorm ####
var = c(1,0.01)
res2 = profile.Cproc(make.Circle(),data,t,pars=pars,coefs,bbasis,var=var, out.meth='BFGS', in.meth='nlminb',control.in=control.in,control.out=control.out)

### SSE####
res1 = profile.sse(make.Circle(),data,t,pars=spars,coefs,bbasis,lambda=lambda,in.meth='house',
out.meth='nlminb',control.in=control.in,control.out=control.out)

#### EM####
proc$more$more$v.more = list(mat=NULL,sub=matrix(c(1,1,3,2,2,3),2,3,byrow=T))
lik$more$v.more = list(mat=NULL,sub=matrix(c(1,1,4,2,2,4),2,3,byrow=T))

pars2 = c(pars,log(c(1e-4,0.0025)))

trfn = function(pars){ return( c(pars[1:2],exp(pars[3:4])) ) }

res2 = EM(pars2,t,data,ncoefs,lik,proc,50,control,tfn=trfn)



### A little experiment

par2 = -c(0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
sses = 0*par2
for(i in 1:length(par2)){
  spars = pars
  spars[2] = par2[i]
  
  res0 = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
        control=control.out,times=t,data=data,lik=lik,proc=proc,pars=spars)

  ncoefs = matrix(res0$par,102,2)
  devals = lik$bvals%*%ncoefs
  
  sses[i] = sum( (data-devals)^2 )
}



## Let's try

profile.obj = sse.setup(pars=pars,coefs,make.genlin(),basisvals=bbasis,lambda,fd.obj=NULL,
    more=list(mat=matrix(0,2,2),sub=matrix(c(1,2,1,2,1,2),2,3,byrow=TRUE)),
    weights=NULL,times=t,quadrature=NULL,eps=1e-6)
lik = profile.obj$lik
proc = profile.obj$proc
coefs = profile.obj$coefs


res0 = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
        control=control.out,times=t,data=data,lik=lik,proc=proc,pars=spars)

ncoefs = matrix(res0$par,102,2)



res1 = profile.sse(make.genlin(),data,t,pars=spars,coefs,bbasis,lambda=lambda,in.meth='nlminb',
more=list(mat=matrix(0,2,2),sub=matrix(c(1,2,1,2,1,2),2,3,byrow=TRUE)),
out.meth='house',control.in=control.in,control.out=control.out)
