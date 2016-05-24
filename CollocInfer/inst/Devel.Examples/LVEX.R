library('CollocInfer')
source('../../R/ls.shortcut.r')
source('LV.R')

#data  = read.table('LynxHare.dat',head=TRUE)
data = read.table('LynxHare2.dat',head=TRUE)
times = data[,1]
data = as.matrix(data[,2:3])



plot(data,type='l')

matplot(times,data,type='b')

bbasis = create.bspline.basis(range(times),norder=4,nbasis=40)
fdobj = smooth.basis(times,data,fdPar(bbasis,2,lambda=0.1))

plotfit.fd(data,times,fdobj$fd)


coefs = fdobj$fd$coefs

fn = make.LV()

setup.obj = LS.setup(rep(0,4),fn=fn,basisvals=bbasis,coefs=coefs,lambda=1e5,data=data,times=times)
lik = setup.obj$lik
proc = setup.obj$proc

proc$more$weights = 1e5 + 0*proc$more$weights


# First some parameters

fvals = eval.fd(proc$more$qpts,fdobj$fd)
Y = eval.fd(proc$more$qpts,fdobj$fd,1)

X1 = cbind(fvals[,1],-fvals[,1]*fvals[,2])
p1 = solve(t(X1)%*%X1,t(X1)%*%Y[,1])

X2 = cbind(-fvals[,2],fvals[,1]*fvals[,2])
p2 = solve(t(X2)%*%X2,t(X2)%*%Y[,2])

pars = c(p1,p2)

#resP = ParsMatchOpt(rep(0,4),coefs,proc)
#pars = resP$pars

# Let's look at the trajectories

y0 = as.vector(eval.fd(times[1],fdobj$fd))
out = lsoda(y0,times,oderhs,parms=list(proc=proc,pars=pars))
matplot(out[,1],out[,2:3],type='l',lty=1)
matplot(times,data,add=TRUE,pch=1,type='b',lty=2)

# Now let's try some trajectory matching

objectivefn = function(pars)
{
  print(pars)
  parms = list(proc=proc,pars=pars[3:6])
  out = lsoda(pars[1:2],times,oderhs,parms)
  traj = out[,2:3]
  val = sum(lik$fn(data,times,traj,pars,lik$more))
  print(val)
  return(val)
}

resT = optim(c(y0,pars),objectivefn,lower=0,method='L-BFGS-B')

y0T = resT$par[1:2]
parsT = resT$par[3:6]

out = lsoda(y0T,proc$more$qpts,oderhs,parms=list(proc=proc,pars=parsT))
matplot(out[,1],out[,2:3],type='l',lty=1)
matplot(times,data,add=TRUE,pch=1,type='b',lty=2)


# Now some profiling

resI = inneropt(data,times,pars,coefs,lik,proc)

matplot(proc$more$qpts,proc$bvals$bvals%*%resI$coefs,type='l',lty=1)
matplot(times,data,type='b',pch=1,lty=2,add=TRUE)


resO = outeropt(data,times,pars,resI$coefs,lik,proc,out.meth='optim')

matplot(times,lik$bvals%*%resO$coefs,type='l',lty=1)
matplot(times,data,type='b',pch=1,lty=2,add=TRUE)


fd5 = fd(resO$coefs,bbasis)
plotfit.fd(data,times,fd5)


### and up the weight factor

proc2 = proc
proc2$more$weights = 1e12 + proc$more$weights*0

resI2 = inneropt(data,times,pars,coefs,lik,proc2)

matplot(proc$more$qpts,proc$bvals$bvals%*%resI2$coefs,type='l',lty=1)
matplot(times,data,type='b',pch=1,lty=2,add=TRUE)

resO2 = outeropt(data,times,pars,resI2$coefs,lik,proc2,out.meth='optim')

matplot(proc$more$qpts,proc$bvals$bvals%*%resO2$coefs,type='l',lty=1)
matplot(times,data,type='b',pch=1,lty=2,add=TRUE)

fd12 = fd(resO2$coefs,bbasis)
plotfit.fd(data,times,fd12)


### Now lets try a larger basis

bbasis = create.bspline.basis(range(times),norder=4,nbasis=200)

setup.obj = LS.setup(rep(0,4),fn=fn,basisvals=bbasis,coefs=coefs,lambda=1e5,data=data,times=times)
lik = setup.obj$lik
proc = setup.obj$proc


fdobj = smooth.basis(proc$more$qpts,eval.fd(proc$more$qpts,fd5),fdPar(bbasis,2,0.001))

coefs = fdobj$fd$coefs

resO.2 = outeropt(data,times,resO$pars,coefs,lik,proc,out.meth='optim')

matplot(proc$more$qpts,proc$bvals$bvals%*%resO.2$coefs,type='l',lty=1)
matplot(times,data,type='b',pch=1,lty=2,add=TRUE)


fdobj = smooth.basis(proc$more$qpts,eval.fd(proc$more$qpts,fd12),fdPar(bbasis,2,0.001))
proc2 = proc
proc2$more$weights = 1e12/5 + proc$more$weights*0


coefs = fdobj$fd$coefs

res2O.2 = outeropt(data,times,resO2$pars,coefs,lik,proc2,out.meth='optim')

matplot(proc$more$qpts,proc$bvals$bvals%*%res2O.2$coefs,type='l',lty=1)
matplot(times,data,type='b',pch=1,lty=2,add=TRUE)
