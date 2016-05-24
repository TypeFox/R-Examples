#### Henon map examples

library('fda')
library('maxLik')
library('MASS')
library('Matrix')
library('SparseM')
source('Henon.R')
source('ProfileR.R')
source('Dproc.R')
source('Multinorm.R')
source('genlin.R')
source('cvar.R')
source('makeid.R')
source('SSEproc.R')
source('SSElik.R')
source('inneropt.R')


hpars = c(1.4,0.3)

ntimes = 200

x = c(-1,1)
X = matrix(0,ntimes+20,2)
X[1,] = x

for(i in 2:(ntimes+20)){ X[i,] = make.Henon()$ode(i,X[i-1,],hpars,NULL) }

X = X[20+1:ntimes,]

Y = X + 0.05*matrix(rnorm(ntimes*2),ntimes,2)

t = 1:ntimes

# Define a discrete-time basis

bbasis = create.bspline.basis(c(0.5,200.5),nbasis=200,norder=1,breaks=seq(0.5,200.5,1))
basisvals = Matrix(eval.basis(1:ntimes,bbasis),sparse=TRUE)


# Now lets define a process likelihood:

proc = make.Dproc()                            # Use discrete-time process likelihood
proc$bvals = basisvals                         # I need basis values 

proc$more = make.multinorm()                   # x_{t+1}|F(x_{t}) is multivariate normalproc$more$qpts = t[1:(ntimes-1)]
                                               # Evaluation times
proc$more$more = c(make.Henon(),make.cvar())   # F = Henon map, variance is constant
proc$more$more$f.more = NULL                   # I don't need any extra inputs into F
proc$more$more$v.more = list(mat=0.04*diag(rep(1,2)),sub=matrix(0,0,3))
                                               # covariance is known diagonal
proc$more$qpts = 1:199


# And an observation likelihood:

lik = make.multinorm()
lik$bvals = basisvals
lik$more = c(make.id(),make.cvar())
lik$more$f.more = NULL
lik$more$v.more = list(mat=diag(rep(1,2)),sub=matrix(0,0,3))


# I need some optimization control variables

control=list(trace = 0,maxit = 1000,maxtry = 10,reltol = 1e-6,meth = "BFGS")

control.in = control
control.in$reltol = 1e-12

control.out = control
control.out$trace = 2


# First thing is to generate a smooth of the data:

coefs = as.vector(Y)

res = optim(coefs,SplineCoefsErr,gr=SplineCoefsDC,hessian=T,   
        method="BFGS",control=control.out,
        times=t,data=Y,lik=lik,proc=proc,pars=hpars)

res1 = SplineEst.NewtRaph(coefs,t,data=Y,lik,proc,hpars,control=control.out)


res2 = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
        control=control.out,times=t,data=Y,lik=lik,proc=proc,pars=hpars)

ncoefs = matrix(res2$par,200,2)

#ProfileEnv = new.env()
#assign('optcoefs',ncoefs,3,ProfileEnv)
#assign('curcoefs',ncoefs,3,ProfileEnv)

f = ProfileErr(hpars,hpars,t,X,ncoefs,lik,proc,in.meth="nlminb",control.in=control)  # outer objective
g = ProfileDP(hpars,hpars,t,X,ncoefs,lik,proc)                             # outer gradient

# Comparison to SSE

sproc = make.SSEproc()
sproc$bvals = list(bvals=basisvals[1:(ntimes-1),],dbvals=basisvals[2:ntimes,])
sproc$more = make.Henon()
sproc$more$weights = matrix(1e6,ntimes-1,2)
sproc$more$qpts = 1:(ntimes-1)

slik = make.SSElik()
slik$more = make.id()
slik$more$weights = matrix(1,ntimes,2)
slik$bvals = basisvals

res2 = SplineEst.NewtRaph(coefs,t,Y,slik,sproc,hpars,control=control.out)


res3 = optim(coefs,SplineCoefsErr,gr=SplineCoefsDC,hessian=T,   
          method="BFGS",control=control.out,
          times=t,data=Y,lik=slik,proc=sproc,pars=hpars)

res4 = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
        control=control.out,times=t,data=Y,lik=slik,proc=sproc,pars=hpars)


# Now we'll estimate some parameters

res5 = nlminb(hpars+c(0.4,0.1),ProfileErr,allpars=hpars+c(0.4,0.1),times=t,data=Y,coefs=ncoefs,lik=lik,proc=proc,
    control.in=control.in,control=control.out,gr=ProfileDP,in.meth="nlminb")



# As an alternative, consider squared error

res6 = Profile.GausNewt(hpars+c(0.4,0.1),t,Y,coefs=ncoefs,slik,sproc,in.meth="nlminb",control.in=control.in)

res7 = optim(hpars+c(0.4,0.1),ProfileErr,allpars=hpars+c(0.4,0.1),times=t,data=Y,coefs=ncoefs,lik=slik,proc=sproc,
    control.in=control.in,control=control.in,gr=ProfileDP,method="BFGS")

res8 = nlminb(hpars,ProfileErr,allpars=hpars+c(0.4,0.1),times=t,data=Y,coefs=ncoefs,lik=slik,proc=sproc,
    control.in=control.in,control=control.out,gradient=ProfileDP)


res9 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
    data = list(times=t,data=Y,coefs=ncoefs,lik=slik,proc=proc,
    in.meth='nlminb',control.in=control.in),
    start = list(pars=hpars),trace=TRUE,control=control.out)

g = res9$m$gradient()
df = diag(res9$m$resid())%*%g
C6 = NeweyWest.Var(t(g)%*%g, df[1:200,]+df[201:400,],5)




# lets make the system stochastic

x = c(-1,1)
X = matrix(0,ntimes+20,2)
X[1,] = x

for(i in 2:(ntimes+20)){ X[i,] = Henon.ode(i,X[i-1,],hpars,NULL) + 0.01*rnorm(2) }

X = X[20+1:ntimes,]

Y = X + 0.05*matrix(rnorm(ntimes*2),ntimes,2)

proc$more$more$v.more$mat = diag(rep(0.001,2))

res10 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
    data = list(times=t,data=Y,coefs=ncoefs,lik=slik,proc=proc,
    in.meth='nlminb',control.in=control.in),
    start = list(pars=hpars),trace=TRUE)



res11 = nlminb(hpars,ProfileErr,allpars=hpars,times=t,data=Y,coefs=ncoefs,lik=slik,proc=sproc,
    control.in=control.in,control=control.out,gradient=ProfileDP)


g = res10$m$gradient()
df = diag(res10$m$resid())%*%g
C10 = NeweyWest.Var(t(g)%*%g, df[1:200,]+df[201:400,],5)



# Now lets try EM

proc$more$more$v.more = list(mat = diag(rep(1e-4,2)))
lik$more$v.more = list(mat = diag(rep(0.0025,2)))

res12 = nlminb(ncoefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
        control=control,times=t,data=Y,lik=lik,proc=proc,pars=hpars)

res13 = optim(hpars,Qfn,method="BFGS",control=control,times=t,data=Y,coefs=res$par,
    lik=lik,proc=proc,oldpars=hpars)


proc$more$more$v.more = list(mat=NULL,sub=matrix(c(1,1,3,2,2,3),2,3,byrow=T))
lik$more$v.more = list(mat=NULL,sub=matrix(c(1,1,4,2,2,4),2,3,byrow=T))

hpars2 = c(hpars,log(c(1e-4,0.0025)))

trfn = function(pars){ return( c(pars[1:2],exp(pars[3:4])) ) }

res14 = EM(hpars2,t,Y,ncoefs,lik,proc,50,control,tfn=trfn)

