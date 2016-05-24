pomp.skeleton = function(times,y,p,more)
{
  # Turns a skeleton function from a 'pomp' object into the right hand side
  # of and ODE for use in CollocInfer
  x <- array(t(y),dim=c(ncol(y),1,nrow(y)),dimnames=list(colnames(y),NULL,NULL))
  params <- array(data=p,dim=c(length(p),1),dimnames=list(names(p),NULL))
  
  y = skeleton(more$pomp.obj,x=x,params=params,t=times)
  return(t(y[,1,]))
}




###

data(ricker)
x <- array(
           data=states(ricker),
           dim=c(2,3,52),
           dimnames=list(rownames(states(ricker)),NULL,NULL)
           )
p <- array(
           data=coef(ricker),
           dim=c(5,3),
           dimnames=list(names(coef(ricker)),NULL)
           )
p["log.r",] <- c(1,2,4)


x = trajectory(ricker)
x = t(x[,1])
p = p[,1]

times=time(ricker,t0=T)

bvals = create.bspline.basis(rangeval=range(times),norder=2,nbasis=length(times))

objs = LS.setup(p,x,pomp.skeleton,100,discrete=1,times = times,basisvals=bvals)

proc = objs$proc
lik = objs$lik

proc$more$more$more = list(pomp.obj=ricker)


res1 = inneropt(x,times,p,x,lik,proc,in.meth='maxNR',control.in=list(trace=2))


#### SIR


p <- array(
           data=coef(euler.sir),
           dim=c(15,3),
           dimnames=list(names(coef(euler.sir)),NULL)
           )
p["beta2",1:2] <- c(3,5)
params=p[,3,drop=TRUE]

x <- trajectory(euler.sir)
x = t(x[,1,])

Y = x + 0.5*matrix(rnorm(length(x)),nrow(x),ncol(x))

times=time(euler.sir,t0=TRUE)

## Still need a basis

bbasis = create.bspline.basis(range(times),norder=4,breaks=seq(0,4,len=51))
mids = seq(0.04,3.96,by=0.08)

obvals = eval.basis(times,bbasis)
bvals = list(bvals = eval.basis(mids,bbasis),dbvals = eval.basis(mids,bbasis,1))


lik = make.SSElik()
lik$more = make.id()
lik$more$weights = matrix(1,nrow(Y),ncol(Y))
lik$bvals = obvals


proc = make.SSEproc()
proc$bvals= bvals
proc$more$weights = 100*matrix(1,nrow(Y),ncol(Y))
proc$more = make.findif.ode()
proc$more$more = list(eps = 1e-6,fn=pomp.skeleton)
proc$more$more$more = list(pomp.obj = euler.sir)


objs = LS.setup(params,coefs,euler.sir,100,discrete=0,times = times,basisvals=bbasis)

proc = objs$proc
lik = objs$lik
proc$more$more$more = list(pomp.obj = euler.sir)



coefs = solve(t(obvals)%*%obvals,t(obvals)%*%Y)

res1 = inneropt(Y,times,params,coefs,lik,proc,in.meth='nlminb')

res2 = outeropt(Y,times,params,coefs,lik,proc,in.meth='nlminb',out.meth='nlminb')

