estdiweibull <-
function(x, method="P", control=list())
{
par<-numeric(2)
n<-length(x)
m1<-mean(x)
m2<-mean(x^2)
beta0<-4
q0<-(m1-1)/m1

# method of proportion
if (method=="P")
{
y<-sum(x==1)

if(y==0)
{
message("Method of proportion not applicable for estimating q and beta!")
par<-c(NA,NA)
}

else
{
par[1]<-y/n
z<-sum(x<=2)

if(z==y | z==round(n))
{
message("Method of proportion not applicable for estimating beta!")
par[2]<-NA
}
else
par[2]<--1/log(2)*log(log(z/n)/log(y/n))
}
}

# method of moments
if(method=="M")
{
if(sum(x<=2)==n)
{
message("Method of moments not applicable on this sample!")
par<-c(NA,NA)
}
else
{
eps<-ifelse(is.null(control$eps),0.0001,control$eps)
nmax<-ifelse(is.null(control$nmax),1000,control$nmax)
par<-solnp(pars<-c(q0,beta0),fun=lossdiw,x=x,eps=eps,nmax=nmax,LB=c(0,2),UB=c(1,100))$pars
}
}

# heuristic/ML
if (method=="H")
{
if(sum(x<=2)==n)
{
message("Heuristic method not applicable on this sample!")
par<-c(NA,NA)
}
else
{
beta1<-ifelse(is.null(control$beta1),1,control$beta1)
z<-ifelse(is.null(control$z),0.1,control$z)
r<-ifelse(is.null(control$r),0.1,control$r)
Leps<-ifelse(is.null(control$Leps),0.0001,control$Leps)
par<-heuristic(x, beta1=beta1, z=z, r=r, Leps=Leps)
}
}

# paper plot (graphical method)
else if (method=="PP")
{
x1<-log(x)
F<-ecdf(x)
y1<--log(-log((F(x)-0.5/n)))
mod<-lm(y1~x1)
beta<-mod$coefficients[[2]]
q<-exp(-exp(-mod$coefficients[[1]]))
par<-c(q,beta)
}
par
}

