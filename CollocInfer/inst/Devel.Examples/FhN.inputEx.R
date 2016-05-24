FhNtimes = seq(0,1000,25)
FhNpars =  c(0.7,0.8,0.008)
names(FhNpars) = c('a','b','c')

knots = seq(0,1000,2)
norder = 3
nbasis = length(knots) + norder - 2

bbasis = create.bspline.basis(range=range(FhNtimes),nbasis=nbasis,
	norder=norder,breaks=knots)
	
coefs = matrix(0,bbasis$nbasis,2)
colnames(coefs) = c('V','R')

profile.obj = LS.setup(pars=FhNpars,coefs=coefs,fn=make.fhn.input()$fn,basisvals=bbasis,lambda=1e8,times=FhNtimes,
                  more = list(infn=function(t){return(0.2*(t>500))}))
lik = profile.obj$lik
proc= profile.obj$proc

# Let's start off with a smooth

FhNdata = DE2x(c(-1,1),FhNtimes,FhNpars,proc)
matplot(FhNtimes,FhNdata,type='b')

fd.obj = smooth.basis(FhNtimes,FhNdata,fdPar(bbasis,1,0.01))
plotfit.fd(FhNdata,FhNtimes,fd.obj$fd)

coefs = fd.obj$fd$coefs

# Now we can look at doing the inner optimization

resI = inneropt(FhNdata,FhNtimes,FhNpars,coefs,lik,proc,in.meth="SplineEst")

traj = proc$bvals$bvals%*%resI$coefs
matplot(proc$more$qpts,traj,type='l')
matplot(FhNtimes,FhNdata,add=TRUE)

# And the outer optimization

resO =  outeropt(FhNdata,FhNtimes,FhNpars,coefs,lik,proc)

traj = proc$bvals$bvals%*%resO$coefs
matplot(proc$more$qpts,traj,type='l')
matplot(FhNtimes,FhNdata,add=TRUE)

dtraj = proc$bvals$dbvals%*%resO$coefs
ftraj = proc$more$fn(proc$more$qpts,traj,FhNpars,proc$more$more)

matplot(proc$more$qpts,dtraj,type='l',lty=1)
matplot(proc$more$qpts,ftraj,type='l',lty=2,add=TRUE)