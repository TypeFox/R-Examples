library('CollocInfer')

# The per-individual growth rates for the Rosenzweig-Macarthur model are given below:

RosMac = function(t,x,p,more)
{
  p = exp(p); x = exp(x)
  dx = x;
  
  dx[,'C'] = p['rho']*(1- x[,'C']/p['kappaC']) - p['gamma']*x[,'B']/(p['kappaB']+x[,'C']) 
  dx[,'B'] = p['chi']*p['gamma']*x[,'C']/(p['kappaB']+x[,'C']) - p['delta']

  return(dx)
}


# We can now obtain some data (from Becks et. al., 2010, Ecology Letters)

data(ChemoRMData)

time = ChemoRMTime
data = log(ChemoRMData)

par(mar=c(5,5,1,1))
matplot(time,data,pch=c('C','B'),xlab='days',ylab='log(Obs)',cex.lab=2.5,cex.axis=2.5,cex=1.5)


# and set up some names for the variables and the parameters

varnames = RMvarnames
parnames = RMparnames

# Define some asis functions

rr = range(time)
knots = seq(rr[1],rr[2],by=1)

mids = c(min(knots),(knots[1:(length(knots)-1)] + 0.5),max(knots))

bbasis = create.bspline.basis(rr,norder=4,breaks=knots)

# An inital set of parameters and coefficients

coef0 = smooth.basis(time,data,fdPar(bbasis,int2Lfd(2),10))$fd$coef
colnames(coef0) = varnames

#pars = c(log(c(0.4,5e5)),rep(0,4))
pars = log(ChemoRMPars)
names(pars) = parnames

# Whic parameters to we wish to estimate?

activepars = 3:6


## Get profiling objects

out = LS.setup(pars=pars,coefs=coef0,basisvals=bbasis,fn=RosMac,lambda=c(5e4,5e2),
          times=time)
          
          
lik = out$lik
proc = out$proc


## Gradient Matching to Start With

res1 = ParsMatchOpt(pars,coef0,proc,active=3:6)


## And profile

res2 = outeropt(data,time,res1$pars,coef0,lik,proc,active=activepars)

# plot the results

out2 = CollocInferPlots(res2$coefs,res2$pars,lik,proc,times=time,data=data,
          cex.lab=2.5,cex.axis=2.5,cex=1.5,lwd=3)

## Now we can look at the covariance

covar = Profile.covariance(res2$pars,times=time,data=data,
          coefs=res2$coefs,lik=lik,proc=proc,active=activepars)
CIs = cbind( res2$pars[activepars] - 2*sqrt(diag(covar)),
               res2$pars[activepars] + 2*sqrt(diag(covar)) )
rownames(CIs) = parnames[activepars]
exp(CIs)

# We can also examine forward prediction error

whichtimes = cbind(1:102,6:107)

lambdafac = c(0.1,0.5,1,5,10)
FPEs = 0*lambdafac
for(ilam in 1:length(lambdafac)){
  t.res = Profile.LS(RosMac,data,time,res2$pars,res2$coefs,bbasis,
                      lambdafac[ilam]*c(5e4,5e2),active=activepars,out.meth='nlminb')
  FPEs[ilam] = forward.prediction.error(time,data,t.res$coefs,
                                   lik,proc,t.res$pars,whichtimes)
}
FPEs





