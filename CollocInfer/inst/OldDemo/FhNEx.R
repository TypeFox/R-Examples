# An example file fitting the FitzHugh-Nagumo equations to data in the
# new R Profiling Code. This will eventually be interfaced with the new R 
# "Partially Observed Markov Process (pomp)" class of objects. 

library('fda')
library('odesolve')
library('maxLik')
library('MASS')
library('Matrix')
library('SparseM')
source('ProfileR.R')
source('sse.shortcut.R')
source('fhn.R')
source('findif.ode.R')
source('SSElik.R')
source('SSEproc.R')
source('makeid.R')
source('inneropt.R')

# First Create Some Data

t = seq(0,20,0.05)

pars = c(0.2,0.2,3)
names(pars) = c('a','b','c')

x0 = c(-1,1)
names(x0)= c('V','R')
y = lsoda(x0,times=t,func=make.fhn()$fn.ode,pars)
y = y[,2:3]

data = y + matrix(rnorm(802),401,2)



# Now a basis object

knots = seq(0,20,0.2)
norder = 3
nbasis = length(knots) + norder - 2
range = c(0,20)

bbasis = create.bspline.basis(range=range,nbasis=nbasis,
	norder=norder,breaks=knots)


# Initial values for coefficients will be obtained by smoothing

fd.data = array(data,c(dim(data)[1],1,dim(data)[2]))

DEfd = data2fd( fd.data,t,bbasis,fdnames=list(NULL,NULL,c('V','R')) )

coefs = matrix(DEfd$coefs,dim(DEfd$coefs)[1],dim(DEfd$coefs)[3])
colnames(coefs) = DEfd$fdnames[[3]]

# Usual meta-parameters; quadrature points, weights and knots

lambda = c(10000,10000)
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
#lik$bvals = as.matrix.csr(eval.basis(t,bbasis))
lik$bvals = Matrix(eval.basis(t,bbasis),sparse=TRUE)

# Proc is a process log likelihood -- in this case treated as squared
# discrepancy from the ODE definition. 

procmore = make.fhn()
procmore$names = varnames
procmore$parnames = parnames

procmore$weights = qwts
procmore$qpts = qpts


proc = make.SSEproc()
proc$more = procmore
#proc$bvals = list(bvals=as.matrix.csr(eval.basis(procmore$qpts,bbasis,0)),
#		dbvals = as.matrix.csr(eval.basis(procmore$qpts,bbasis,1)))
proc$bvals = list(bvals=Matrix(eval.basis(procmore$qpts,bbasis,0),sparse=TRUE),
		dbvals = Matrix(eval.basis(procmore$qpts,bbasis,1),sparse=TRUE))




###Now lets try some optimization

spars = c(0.2,0.2,2)          # Perturbed parameters


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

# We'll try a simple SSE setup:


res = sse.setup(pars=pars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,lambda=lambda,times=t)

res = smooth.sse(make.fhn(),data,t,pars,coefs,bbasis,lambda=lambda,control.in=control.in)

res = profile.sse(make.fhn(),data,t,pars,coefs,bbasis,lambda=lambda,out.meth='nls',
	control.in=control.in,control.out=control.out)


# Alternative is simply to use the functional data object


res = sse.setup(pars=pars,fd.obj=DEfd,fn=make.fhn(),lambda=lambda,times=t)

res = smooth.sse(pars=pars,fd.obj=DEfd,fn=make.fhn(),lambda=lambda,times=t,data=fd.data,control.in=control.in)

res = profile.sse(fn=make.fhn(),data=fd.data,times=t,pars=pars,fd.obj=DEfd,lambda=lambda,out.meth='house',
	control.in=control.in,control.out=control.out)


# Start with the inner optimization, with various optimization methods

f = SplineCoefsErr(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)
g = SplineCoefsDC(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)
h = SplineCoefsDC2(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)
g2 = SplineCoefsDP(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)
h2 = SplineCoefsDCDP(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)


res = optim(coefs,SplineCoefsErr,gr=SplineCoefsDC,hessian=T,   		
		method="BFGS",control=control.out,
		times=t,data=data,lik=lik,proc=proc,pars=spars)

res0 = nlminb(coefs+2,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
	control=control.out,times=t,data=data,lik=lik,proc=proc,pars=spars)


res2 = maxNR(SplineCoefsErr,start=as.vector(coefs),times=t,data=data,lik=lik,proc=proc,pars=spars,sgn=-1,
	grad=SplineCoefsDC,hess=SplineCoefsDC2,print.level=2,iterlim=100)


res3 = SplineEst.NewtRaph(coefs,t,data,lik,proc,spars)


# record the coefficients for the sake of good starting values

ncoefs = array(res0$par,c(bbasis$nbasis,ncol(data)))

ProfileEnv = new.env()
assign('optcoef',ncoefs,3,ProfileEnv)
assign('curcoefs',ncoefs,3,ProfileEnv)

# Test outer objective criterion

f = ProfileErr(pars,pars,t,data,ncoefs,lik,proc,in.meth="nlminb",control.in=control.in)  
g = ProfileDP(pars,pars,t,data,ncoefs,lik,proc,sum=FALSE)                            


# There are specific criteria for SSE so that a Gauss-Newton iteration can be used

res4 = ProfileSSE(pars,t,data,ncoefs,lik,proc,in.meth="nlminb",control.in) 

# Gauss-Newton function coded in R. 

res5 = Profile.GausNewt(spars,t,data,ncoefs,lik,proc,in.meth="nlminb",control.in)
 
# the alternative standard function is NLS

res6 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth="nlminb",control.in=control.in),start = list(pars=pars),
	trace=TRUE)

# for which we can try a Newey-West type covariance estimate

g = res6$m$gradient()
df = diag(res6$m$resid())%*%g
C6 = NeweyWest.Var(t(g)%*%g, df[1:401,]+df[402:802,],5)


# Try using a simple second-derivative optimizer for the inner optimization:

res7 = Profile.GausNewt(spars,t,data,ncoefs,lik,proc,in.meth='house',control.in)

# nls with a number of alternative inner optimizations: self defined:
 
res8 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='house',control.in=control.in),
	start = list(pars=spars),trace=TRUE)

# nlminb:


res9 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='nlminb',control.in=control.in),
	start = list(pars=spars),trace=TRUE)

# maxNR:


res10 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='maxNR',control.in=control.in),
	start = list(pars=spars),trace=TRUE)



# When squared error is not being employed, things proceed somewhat more slowly
# However, this allows the smoothing criteria to fit into the pomp framework
# fairly readily. 


res11 = optim(pars,ProfileErr,allpars=pars,times=t,data=data,coef=ncoefs,lik=lik,proc=proc,hessian=T,
	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP,method="BFGS")

g = ProfileDP(res11$par,res11$par,t,data,ncoefs,lik,proc,sum=FALSE)  

gg = apply(g,2,sum)

H = 0*res11$hess

for(i in 1:length(pars)){
	tpars = res11$par
	tpars[i] = tpars[i] + 1e-4

	tf = ProfileErr(tpars,tpars,t,data,ncoefs,lik,proc,in.meth="nlminb",control.in=control.in)  
	tg = ProfileDP(tpars,tpars,t,data,ncoefs,lik,proc,sum=TRUE)                            

	H[,i] = (tg-gg)*1e4
}
	

C11 = NeweyWest.Var(H,g,5)


# optim searches a much larger space, which can be a problem, but it finds a better
# minimum than 

res12 = nlminb(pars,ProfileErr,allpars=pars,times=t,data=data,coef=ncoefs,lik=lik,proc=proc,
	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP)


# Now lets see what happens for a finite-difference approximation

proc2more = make.findif.ode()
proc2more$more = list()
proc2more$more$fn = fhn.fun
proc2more$names = varnames
proc2more$parnames = parnames

proc2more$more$eps = 1e-6

proc2more$weights = qwts
proc2more$qpts = qpts

proc2 = make.SSEproc()
proc2$more = proc2more
proc2$bvals = list(bvals=eval.basis(procmore$qpts,bbasis,0),
		dbvals = eval.basis(procmore$qpts,bbasis,1))


res13 = Profile.GausNewt(pars,t,data,ncoefs,lik,proc2,in.meth='nlminb',control.in)

res14 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc2,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='nlminb',control.in=control.in),
	start = list(pars=spars),trace=TRUE)



# Now we'll only observe the first component

data2 = data
data2[,2] = NA

res15 = Profile.GausNewt(pars,t,data2,ncoefs,lik,proc,in.meth='nlminb',control.in)

res16 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc2,in.meth,control.in),
	data = list(times=t,data=data2,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='nlminb',control.in=control.in),
	start = list(pars=spars),trace=TRUE)



# Lets suppose that I've measured a linear combination

source('genlin.R')

lik2more = make.genlin()
lik2more$more = list()
lik2more$weights = weights[,1]

lik2more$more = list()
lik2more$more$mat = matrix(c(1,0),1,2)
lik2more$more$sub = matrix(0,0,3)

lik2 = lik;
lik2$more = lik2more;

data3 = as.matrix(data[,1],nrow(data),1)

res17 = Profile.GausNewt(pars,t,data3,ncoefs,lik2,proc,in.meth='house',control.in)

# Now lets add estimating a set of parameters in lik to the mix

lik3more = lik2more
lik3more$more$mat = matrix(0,1,2)
lik3more$more$sub = matrix(c(1,1,4,1,2,5),2,3,byrow=T)

lik3 = lik2
lik3$more = lik3more

proc2 = proc
proc2$more$parnames = c(proc$more$parnames,'l1','l2')

pars2 = c(pars,1,0)
names(pars2) = proc2$more$parnames

res18 = Profile.GausNewt(pars2,t,data3,ncoefs,lik3,proc2,in.meth='nlminb',control.in)


res19 = optim(pars2,ProfileErr,allpars=pars2,times=t,data=data3,ncoefs,lik=lik3,proc=proc2,
	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP,method="BFGS")


res20 = nlminb(pars2,ProfileErr,allpars=pars2,times=t,data=data3,coef=ncoefs,lik=lik3,proc=proc2,
	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP)

f = ProfileErr(res20$par,res20$par,t,data3,ncoefs,lik3,proc2,in.meth="nlminb",control.in=control.in)  
g = ProfileDP(res20$par,res20$par,t,data3,ncoefs,lik3,proc2,sum=FALSE)  

gg = apply(g,2,sum)

H = matrix(0,5,5)

for(i in 1:5){
	tpars = res20$par
	tpars[i] = tpars[i] + 1e-4

	tf = ProfileErr(tpars,tpars,t,data3,ncoefs,lik3,proc2,in.meth="nlminb",control.in=control.in)  
	tg = ProfileDP(tpars,tpars,t,data3,ncoefs,lik3,proc2,sum=TRUE)                            

	H[,i] = (tg-gg)*1e4
}

C20 = NeweyWest.Var(0.5*(H+t(H)),g,5)



res21 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data3,coefs=ncoefs,lik=lik3,proc=proc2,
	in.meth='nlminb',control.in=control.in),
	start = list(pars=pars2),trace=TRUE)

ff =  ProfileSSE(res20$par,t,data3,ncoefs,lik3,proc2,in.meth="nlminb",control.in) 


res22 = maxNR(ProfileErr,start=pars2,allpars=pars2,times=t,data=data3,coef=ncoefs,lik=lik3,proc=proc2,
	in.meth='nlminb',control.in=control.in,sgn=-1,grad=ProfileDP,print.level=2,iterlim=100)



# Now we'll set up some repeat experiments

x02 = c(1,-1)
names(x02)= c('V','R')
y2 = lsoda(x02,times=t,func=fhn.fun.ode,pars)
y2 = y2[,2:3]

data2 = y2 + matrix(rnorm(802),401,2)

replik = lik
replik$bvals = diag(rep(1,2))%x%lik$bvals
replik$more$weights = rbind(lik$more$weights,lik$more$weights)

repproc = proc
repproc$bvals = list(bvals = diag(rep(1,2))%x%proc$bvals$bvals,
  dbvals=diag(rep(1,2))%x%proc$bvals$dbvals)
repproc$more$weights = rbind(proc$more$weights,proc$more$weights)

reptimes = c(t,t+max(t))

repdata = rbind(data,data2)

coefs3 = solve( t(replik$bvals)%*%replik$bvals )%*%( t(replik$bvals)%*%repdata )

res23 = nlminb(coefs3,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
	control=control.out,times=reptimes,data=repdata,lik=replik,proc=repproc,pars=spars)


ncoefs = array(res23$par,dim(coefs3))

ProfileEnv = new.env()
assign('optcoef',ncoefs,3,ProfileEnv)
assign('curcoefs',ncoefs,3,ProfileEnv)

res24 = Profile.GausNewt(spars,reptimes,repdata,ncoefs,replik,repproc,in.meth="nlminb",control.in)
 


# Lets try this assuming we have functional data

fd.data2 = array(0,c(nrow(data2),2,2))

fd.data2[,2,] = data2
fd.data2[,1,] = data


DEfd2 = data2fd(fd.data2,t,bbasis,fdnames=list(NULL,NULL,c('V','R')) )

res = sse.setup(pars=pars,fd.obj=DEfd2,fn=make.fhn(),lambda=100,times=t)

res = smooth.sse(pars=pars,fd.obj=DEfd2,fn=make.fhn(),lambda=100,times=t,data=fd.data2,control.in=control.in)

res = profile.sse(fn=make.fhn(),data=fd.data2,times=t,pars=pars,fd.obj=DEfd2,,lambda=100,out.meth='nls',
	control.in=control.in,control.out=control.out)


# Alternatively, we can just do the sse setup thing

coefs2 = DEfd2$coefs


dimnames(coefs2) = list(NULL,NULL,c('V','R'))


res = sse.setup(pars=pars,coefs=coefs,fn=make.fhn(),basisvals=bbasis,lambda=100,times=t)

res = smooth.sse(make.fhn(),data,t,pars,coefs,bbasis,lambda=100,control.in=control.in)

res = profile.sse(make.fhn(),data,t,pars,coefs,bbasis,lambda=100,out.meth='nls',
	control.in=control.in,control.out=control.out)


# Alternative is simply to use the functional data object


res = sse.setup(pars=pars,fd.obj=DEfd,fn=make.fhn(),lambda=100,times=t)

res = smooth.sse(pars=pars,fd.obj=DEfd,fn=make.fhn(),lambda=100,times=t,data=fd.data,control.in=control.in)

res = profile.sse(fn=make.fhn(),data=fd.data,times=t,pars=pars,fd.obj=DEfd,,lambda=100,out.meth='nls',
	control.in=control.in,control.out=control.out)


# Start with the inner optimization, with various optimization methods

f = SplineCoefsErr(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)
g = SplineCoefsDC(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)
h = SplineCoefsDC2(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)
g2 = SplineCoefsDP(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)
h2 = SplineCoefsDCDP(coefs,times=t,data=data,lik=lik,proc=proc,pars=spars)


res = optim(coefs,SplineCoefsErr,gr=SplineCoefsDC,hessian=T,   		
		method="BFGS",control=control.out,
		times=t,data=data,lik=lik,proc=proc,pars=spars)

res0 = nlminb(coefs+2,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
	control=control.out,times=t,data=data,lik=lik,proc=proc,pars=spars)


res2 = maxNR(SplineCoefsErr,start=as.vector(coefs),times=t,data=data,lik=lik,proc=proc,pars=spars,sgn=-1,
	grad=SplineCoefsDC,hess=SplineCoefsDC2,print.level=2,iterlim=100)


res3 = SplineEst.NewtRaph(coefs,t,data,lik,proc,spars)


# record the coefficients for the sake of good starting values

ncoefs = array(res0$par,c(bbasis$nbasis,ncol(data)))

ProfileEnv = new.env()
assign('optcoef',ncoefs,3,ProfileEnv)
assign('curcoefs',ncoefs,3,ProfileEnv)

# Test outer objective criterion

f = ProfileErr(pars,pars,t,data,ncoefs,lik,proc,in.meth="nlminb",control.in=control.in)  
g = ProfileDP(pars,pars,t,data,ncoefs,lik,proc,sum=FALSE)                            


# There are specific criteria for SSE so that a Gauss-Newton iteration can be used

res4 = ProfileSSE(pars,t,data,ncoefs,lik,proc,in.meth="nlminb",control.in) 

# Gauss-Newton function coded in R. 

res5 = Profile.GausNewt(spars,t,data,ncoefs,lik,proc,in.meth="nlminb",control.in)
 
# the alternative standard function is NLS

res6 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth="nlminb",control.in=control.in),start = list(pars=pars),
	trace=TRUE)

# for which we can try a Newey-West type covariance estimate

g = res6$m$gradient()
df = diag(res6$m$resid())%*%g
C6 = NeweyWest.Var(t(g)%*%g, df[1:401,]+df[402:802,],5)


# Try using a simple second-derivative optimizer for the inner optimization:

res7 = Profile.GausNewt(spars,t,data,ncoefs,lik,proc,in.meth='house',control.in)

# nls with a number of alternative inner optimizations: self defined:
 
res8 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='house',control.in=control.in),
	start = list(pars=spars),trace=TRUE)

# nlminb:


res9 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='nlminb',control.in=control.in),
	start = list(pars=spars),trace=TRUE)

# maxNR:


res10 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='maxNR',control.in=control.in),
	start = list(pars=spars),trace=TRUE)



# When squared error is not being employed, things proceed somewhat more slowly
# However, this allows the smoothing criteria to fit into the pomp framework
# fairly readily. 


res11 = optim(pars,ProfileErr,allpars=pars,times=t,data=data,coef=ncoefs,lik=lik,proc=proc,hessian=T,
	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP,method="BFGS")

g = ProfileDP(res11$par,res11$par,t,data,ncoefs,lik,proc,sum=FALSE)  

gg = apply(g,2,sum)

H = 0*res11$hess

for(i in 1:length(pars)){
	tpars = res11$par
	tpars[i] = tpars[i] + 1e-4

	tf = ProfileErr(tpars,tpars,t,data,ncoefs,lik,proc,in.meth="nlminb",control.in=control.in)  
	tg = ProfileDP(tpars,tpars,t,data,ncoefs,lik,proc,sum=TRUE)                            

	H[,i] = (tg-gg)*1e4
}
	

C11 = NeweyWest.Var(H,g,5)


# optim searches a much larger space, which can be a problem, but it finds a better
# minimum than 

res12 = nlminb(pars,ProfileErr,allpars=pars,times=t,data=data,coef=ncoefs,lik=lik,proc=proc,
	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP)


# Now lets see what happens for a finite-difference approximation

proc2more = make.findif.ode()
proc2more$more = list()
proc2more$more$fn = fhn.fun
proc2more$names = varnames
proc2more$parnames = parnames

proc2more$more$eps = 1e-6

proc2more$weights = qwts
proc2more$qpts = qpts

proc2 = make.SSEproc()
proc2$more = proc2more
proc2$bvals = list(bvals=eval.basis(procmore$qpts,bbasis,0),
		dbvals = eval.basis(procmore$qpts,bbasis,1))


res13 = Profile.GausNewt(pars,t,data,ncoefs,lik,proc2,in.meth='nlminb',control.in)

res14 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc2,in.meth,control.in),
	data = list(times=t,data=data,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='nlminb',control.in=control.in),
	start = list(pars=spars),trace=TRUE)



# Now we'll only observe the first component

data2 = data
data2[,2] = NA

res15 = Profile.GausNewt(pars,t,data2,ncoefs,lik,proc,in.meth='nlminb',control.in)

res16 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc2,in.meth,control.in),
	data = list(times=t,data=data2,coefs=ncoefs,lik=lik,proc=proc,
	in.meth='nlminb',control.in=control.in),
	start = list(pars=spars),trace=TRUE)



# Lets suppose that I've measured a linear combination

source('genlin.R')

lik2more = make.genlin()
lik2more$more = list()
lik2more$weights = weights[,1]

lik2more$more = list()
lik2more$more$mat = matrix(c(1,0),1,2)
lik2more$more$sub = matrix(0,0,3)

lik2 = lik;
lik2$more = lik2more;

data3 = as.matrix(data[,1],nrow(data),1)

res17 = Profile.GausNewt(pars,t,data3,ncoefs,lik2,proc,in.meth='house',control.in)

# Now lets add estimating a set of parameters in lik to the mix

lik3more = lik2more
lik3more$more$mat = matrix(0,1,2)
lik3more$more$sub = matrix(c(1,1,4,1,2,5),2,3,byrow=T)

lik3 = lik2
lik3$more = lik3more

proc2 = proc
proc2$more$parnames = c(proc$more$parnames,'l1','l2')

pars2 = c(pars,1,0)
names(pars2) = proc2$more$parnames

res18 = Profile.GausNewt(pars2,t,data3,ncoefs,lik3,proc2,in.meth='nlminb',control.in)


res19 = optim(pars2,ProfileErr,allpars=pars2,times=t,data=data3,ncoefs,lik=lik3,proc=proc2,
	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP,method="BFGS")


res20 = nlminb(pars2,ProfileErr,allpars=pars2,times=t,data=data3,coef=ncoefs,lik=lik3,proc=proc2,
	in.meth='nlminb',control.in=control.in,control=control.out,gr=ProfileDP)

f = ProfileErr(res20$par,res20$par,t,data3,ncoefs,lik3,proc2,in.meth="nlminb",control.in=control.in)  
g = ProfileDP(res20$par,res20$par,t,data3,ncoefs,lik3,proc2,sum=FALSE)  

gg = apply(g,2,sum)

H = matrix(0,5,5)

for(i in 1:5){
	tpars = res20$par
	tpars[i] = tpars[i] + 1e-4

	tf = ProfileErr(tpars,tpars,t,data3,ncoefs,lik3,proc2,in.meth="nlminb",control.in=control.in)  
	tg = ProfileDP(tpars,tpars,t,data3,ncoefs,lik3,proc2,sum=TRUE)                            

	H[,i] = (tg-gg)*1e4
}

C20 = NeweyWest.Var(0.5*(H+t(H)),g,5)



res21 = nls(~ProfileSSE(pars,times,data,coefs,lik,proc,in.meth,control.in),
	data = list(times=t,data=data3,coefs=ncoefs,lik=lik3,proc=proc2,
	in.meth='nlminb',control.in=control.in),
	start = list(pars=pars2),trace=TRUE)

ff =  ProfileSSE(res20$par,t,data3,ncoefs,lik3,proc2,in.meth="nlminb",control.in) 


res22 = maxNR(ProfileErr,start=pars2,allpars=pars2,times=t,data=data3,coef=ncoefs,lik=lik3,proc=proc2,
	in.meth='nlminb',control.in=control.in,sgn=-1,grad=ProfileDP,print.level=2,iterlim=100)



# Now we'll set up some repeat experiments

x02 = c(1,-1)
names(x02)= c('V','R')
y2 = lsoda(x02,times=t,func=fhn.fun.ode,pars)
y2 = y2[,2:3]

data2 = y2 + matrix(rnorm(802),401,2)

replik = lik
replik$bvals = diag(rep(1,2))%x%lik$bvals
replik$more$weights = rbind(lik$more$weights,lik$more$weights)

repproc = proc
repproc$bvals = list(bvals = diag(rep(1,2))%x%proc$bvals$bvals,
  dbvals=diag(rep(1,2))%x%proc$bvals$dbvals)
repproc$more$weights = rbind(proc$more$weights,proc$more$weights)

reptimes = c(t,t+max(t))

repdata = rbind(data,data2)

coefs3 = solve( t(replik$bvals)%*%replik$bvals )%*%( t(replik$bvals)%*%repdata )

res23 = nlminb(coefs3,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
	control=control.out,times=reptimes,data=repdata,lik=replik,proc=repproc,pars=spars)


ncoefs = array(res23$par,dim(coefs3))

ProfileEnv = new.env()
assign('optcoef',ncoefs,3,ProfileEnv)
assign('curcoefs',ncoefs,3,ProfileEnv)

res24 = Profile.GausNewt(spars,reptimes,repdata,ncoefs,replik,repproc,in.meth="nlminb",control.in)
 


# Lets try this assuming we have functional data

fd.data2 = array(0,c(nrow(data2),2,2))

fd.data2[,2,] = data2
fd.data2[,1,] = data


DEfd2 = data2fd(fd.data2,t,bbasis,fdnames=list(NULL,NULL,c('V','R')) )

res = sse.setup(pars=pars,fd.obj=DEfd2,fn=make.fhn(),lambda=lambda,times=t)

res = smooth.sse(pars=pars,fd.obj=DEfd2,fn=make.fhn(),lambda=lambda,times=t,data=fd.data2,control.in=control.in,in.meth='house')

res = profile.sse(fn=make.fhn(),data=fd.data2,times=t,pars=pars,fd.obj=DEfd2,,lambda=lambda,out.meth='nls',
  control.in=control.in,control.out=control.out)


# Alternatively, we can just do the sse setup thing

coefs2 = DEfd2$coefs
dimnames(coefs2) = list(NULL,NULL,c('V','R'))


res = sse.setup(pars=pars,coefs=coefs2,fn=make.fhn(),basisvals=bbasis,lambda=10000,times=t)

res = smooth.sse(fn=make.fhn(),data=fd.data2,times=t,pars=pars,coefs=coefs2,basisvals=bbasis,
  lambda=10000,control.in=control.in)

res = profile.sse(fn=make.fhn(),data=fd.data2,times=t,pars=pars,coefs=coefs2,basisvals=bbasis,
  lambda=1e8,out.meth='nls',control.in=control.in,control.out=control.out,in.meth='house')
