linear.quan<-function(x,y,p=0.5)
{
y<-matrix(y,length(y),1)
n<-dim(x)[1]
d<-dim(x)[2]

rho<-function(t){
    t*(p-(t<0))
}

fn<-function(b) { 
    b2<-matrix(b[2:(d+1)],d,1)
    gx<-b[1]+x%*%b2
    ro<-rho(y-gx)
    return(sum(ro)/n)
}

li<-linear(x,y)
par<-c(li$beta0,li$beta1)    # initial value
par.lower<-rep(-1,d)
par.upper<-rep(1,d)

op.method<-"L-BFGS-B"
op<-optim(par=par,fn=fn,method=op.method)
theta<-op$par

beta0<-theta[1]
beta1<-theta[2:(d+1)]

return(list(beta0=beta0,beta1=beta1))
}


