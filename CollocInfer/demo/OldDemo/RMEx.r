library('CollocInfer')


# Two-species Rosensweig-Macarthur Model

RosMac = function(t,x,p,more){

  p = exp(p)
  dx = x

  dx[,'C1'] = p['r1']*x[,'C1']*(1- x[,'C1']/p['Kc1']- x[,'C2']/p['Kc2']) - p['p']*p['G']*x[,'C1']*x[,'B']/(p['KB']+p['p']*x[,'C1']+x[,'C2'])
  dx[,'C2'] = p['r2']*x[,'C2']*(1- x[,'C1']/p['Kc1']- x[,'C2']/p['Kc2']) - p['G']*x[,'C2']*x[,'B']/(p['KB']+p['p']*x[,'C1']+x[,'C2'])
  dx[,'B'] = p['chiB']*p['G']*(p['p']*x[,'C1']+x[,'C2'])*x[,'B']/(p['KB']+p['p']*x[,'C1']+x[,'C2']) - p['delta']*x[,'B']

  return(dx)
}



# Now I need interesting parameters

RMpars = c(0.2,0.025,0.125,2.2e4,1e5,5e6,1,1e9,0.3)
RMParnames = c('p','r1','r2','Kc1','Kc2','G','chiB','KB','delta') 

logpars = log(RMpars)

names(logpars) = RMParnames

# And initial conditions

RMVarnames = c('C1','C2','B')

x0 = c(50,50,2)
names(x0) = RMVarnames

# Suitable for lsoda

RosMacODE = function(t,z,p){
  p = exp(p)
  x = exp(z)
  dx = x

  dx['C1'] = p['r1']*x['C1']*(1- x['C1']/p['Kc1']-x['C2']/p['Kc2']) - p['p']*p['G']*x['C1']*x['B']/(p['KB']+p['p']*x['C1']+x['C2'])
  dx['C2'] = p['r2']*x['C2']*(1- x['C2']/p['Kc2']- x['C1']/p['Kc1']) - p['G']*x['C2']*x['B']/(p['KB']+p['p']*x['C1']+x['C2'])
  dx['B'] = p['chiB']*p['G']*(p['p']*x['C1']+x['C2'])*x['B']/(p['KB']+p['p']*x['C1']+x['C2']) - p['delta']*x['B']

  return(list(dx/x))
}



# Solve the ODE

time = 0:200

res0 = lsoda(log(x0),time,RosMacODE,p = logpars)
matplot(exp(res0[,2:4]),type='l')

data = res0[,2:4] + 0.2*matrix(rnorm(603),201,3)
colnames(data) = RMVarnames

matplot(data,cex.lab=1.5,cex.axis=1.5)
matplot(res0[,2:4],type='l',add=TRUE)


# And the usual setup


rr = range(time)
knots = seq(rr[1],rr[2],by=1)

bbasis = create.bspline.basis(rr,norder=4,breaks=knots)

# And obtain an initial set of parameters and coefficients

coef0 = smooth.basis(time,data,fdPar(bbasis,int2Lfd(2),10))$fd$coef
colnames(coef0) = RMVarnames


## Get profiling objects

out = LS.setup(pars=logpars,coefs=coef0,basisvals=bbasis,fn=RosMac,lambda=1e5,
          times=time,posproc=TRUE)
          
          
lik = out$lik
proc = out$proc


# New parameter estimates

res1 = ParsMatchOpt(logpars,coef0,proc)
res3 = outeropt(data,time,res1$pars,coef0,lik,proc)
exp(res3$pars)
out3 = CollocInferPlots(res3$coefs,res3$pars,lik,proc,times=time,data=data)


## Now if we only observe the sum of C1 and C2

data2 = cbind( log( exp(data[,'C1'])+exp(data[,'C2'])), data[,'B'])
matplot(data2,cex.lab=1.5,cex.axis=1.5)


RMobsfn = function(t,x,p,more)
{
  x = exp(x)
  y = cbind( x[,'C1']+x[,'C2'],x[,'B'])
  return(log(y))
}

out = LS.setup(pars=logpars,coefs=coef0,basisvals=bbasis,fn=RosMac,lambda=1e5,
          times=time,posproc=TRUE,likfn=RMobsfn)
          
          
lik2 = out$lik
proc2 = out$proc

coef02 = coef0
coef02[,1:2] = 0

Fres3 = FitMatchOpt(coef02,1:2,res1$pars,proc2)
res32 = outeropt(data2,time,res1$pars,Fres3$coefs,lik2,proc2)
exp(res32$pars)
out32 = CollocInferPlots(res32$coefs,res32$pars,lik2,proc2,times=time,data=data2)

### Repeated experiments

x03 = c(15,25,4)
names(x03) = RMVarnames

res03 = lsoda(log(x03),time,RosMacODE,p = logpars)

data03 =  res03[,2:4] + 0.2*matrix(rnorm(603),201,3)

alldat = array(0,c(201,2,3))
alldat[,1,] = data
alldat[,2,] = data03

coef3 = smooth.basis(time,data03,fdPar(bbasis,int2Lfd(2),10))$fd$coef

coefs = array(0,c(dim(coef3)[1],2,3))
coefs[,1,] = coef0
coefs[,2,] = coef3


out = LS.setup(pars=logpars,coefs=coefs,basisvals=bbasis,fn=RosMac,lambda=1e5,
          times=time,data=alldat,posproc=TRUE,names=RMVarnames)
          
          
lik3 = out$lik
proc3 = out$proc

res13 = inneropt(data=out$data,times=out$times,pars=res1$pars,coefs=out$coefs,lik=lik3,proc=proc3)
res33 = outeropt(data=out$data,times=out$times,res1$pars,res13$coefs,lik3,proc3)
exp(res33$pars)
out3 = CollocInferPlots(res33$coefs,res33$pars,lik3,proc3,times=out$times,data=out$data)



