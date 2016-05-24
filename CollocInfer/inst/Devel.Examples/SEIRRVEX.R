source('SEIRRV.R')

library('CollocInfer')
library('pomp.devel')

#source('../../R/OuterOptimization.R')
#source('../../R/InnerOptimization.R')
#source('../../R/ProfileCovariance.R')
#source('../../R/SSEproc.R')

# Extract some data

data(pertussis.sim)
truth = pertussis.sim[["SEIRR.small"]]

pars = coef(truth)
pars = pars[1:16]
pars[1:15] = exp(pars[1:15])

varnames = c('S','E','I','R1','R2')
parnames = names(pars)
active = (1:16)[-c(1:2,16)]

times = time(truth)
data = t(data.array(truth))
truestates = t(states(truth,varnames))

data = data[times<=20,]
truestates = truestates[times < 20,]
times = times[times<=20]

data = matrix(data,length(data),1)

# Basis functions

rr = range(times)
knots = seq(rr[1],rr[2],by=1/52)

mids = c(min(knots),(knots[1:(length(knots)-1)] + 1/104),max(knots))

bbasis = create.bspline.basis(rr,norder=4,breaks=knots)

# proc object

lambda = rep(100,5)

proc = make.SSEproc()
proc$more = make.findif.ode()
proc$more$more = list(fn=make.logtrans()$fn,eps=1e-8)
proc$more$more$more = list(fn=SEIRR)

proc$more$names = varnames
proc$more$parnames = parnames

proc$bvals = list(bvals = Matrix(eval.basis(mids,bbasis),sparse=TRUE),
                 dbvals = Matrix(eval.basis(mids,bbasis,1),sparse=TRUE));
#proc$bvals = list(bvals = eval.basis(mids,bbasis),
#                 dbvals = eval.basis(mids,bbasis,1));
proc$more$qpts = mids
proc$more$weights = lambda

# lik object

lik = make.SSElik()
lik$more = make.genlin()
lik$more$more = list(mat=matrix(0,1,5,byrow=TRUE), sub = matrix(c(1,3,11),1,3,byrow=TRUE))

lik$bvals = Matrix(eval.basis(times,bbasis),sparse=TRUE)
#lik$bvals = eval.basis(times,bbasis)
lik$more$weights = matrix(1,length(times),1)

lik2 = make.logstate.lik()
lik2$more = lik
lik2$bvals = lik$bvals

# Now we can try some smoothing
# start with the true states


DEfd = smooth.basis(times,log(truestates),fdPar(bbasis,1,0.0001))
coefs = DEfd$fd$coefs

res0 = inneropt(coefs=coefs,pars=res1$pars,times=times,data=data,lik=lik2,proc=proc,in.meth='SplineEst',control.in=list(trace=6))

traj = lik$bvals%*%res0$coefs
matplot(times,traj,type='l',lty=1)
matplot(times,log(truestates),type='l',lty=2,add=TRUE)


plot(times,exp(traj[,3])*pars['report.prob'],type='l')
points(times,data)

res1 = outeropt(coefs=res1$coefs,pars=res1$pars,times=times,data=data,lik=lik2,proc=proc,in.meth='SplineEst',control.in=list(trace=1))

traj = lik$bvals%*%res1$coefs
matplot(times,traj,type='l',lty=1)
matplot(times,log(truestates),type='l',lty=2,add=TRUE)


plot(times,exp(traj[,3])*res1$pars['report.prob'],type='l')
points(times,data)