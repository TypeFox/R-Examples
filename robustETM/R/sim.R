#.First.lib <- function(lib, pkg) library.dynam("robustETM", pkg, lib)


sim=function(itr, K, cc, i.n, isetting, lambda, distn)
{
set.seed(itr)

logistic.wald1=function(x, y, h1, h2){
   m=length(x); n=length(y)
   z=c(x, y); h1.z = h1(z); h2.z = h2(z)
   case=c(rep(0,m),rep(1,n))
   fit=glm(case~1+h1.z,family="binomial") 
   temp1 = fit$coefficient[-1]
   temp=t(temp1)%*%solve(vcov(fit)[-1,-1])%*%temp1 
   temp=as.numeric(temp)
   return(temp)
}


###########################################################################################
n.settings=array(NA, c(2, 4))
n.settings[,1]=c(20,20)
n.settings[,2]=c(50,50)
n.settings[,3]=c(100, 100) 
n.settings[,4]=c(500, 500) 

normal.settings=rbind(c(0, 1,0,0.4),
c(1,1,1.7,1.5))


beta.settings=rbind(c(99.04095, 93.75379, 29.58310, 33.86302),
c(24.76013, 19.89309,  7.39583,  7.72513))


gamma.settings=rbind(c(6.25, 10.2400, 2.250000, 2.8224000),
c(0.40,  0.3125, 1.111111, 0.9920635))


t.settings=rbind(c(5.960491, 7.488506, 4.968230,  6.982013),
c(7.086485, 9.286753, 3.188841, 26.231343))

nbinom.settings=rbind(round(c(12.2517616, 40.0624963, 1.5942072 , 2.1432008),0),
c(0.6363854,  0.8167938 ,0.1855292, 0.2185168))
###########################################################################################


n=n.settings[,i.n]
nn=n[1]+n[2]

beta.init=c(0,0)

lamda_length=10;alpha_length=10
lamda_ini=seq(0,1,length=lamda_length)

alpha_ini=seq(-2,2,length=alpha_length)

plrt_ini=plrt_temp=mplrt_ini=matrix(0,nrow=lamda_length,ncol=alpha_length)
lrt_ini=rep(0,lamda_length)


####################################
if (distn=="norm"){
dist.isetting=normal.settings
m1=dist.isetting[1,1]; s1=dist.isetting[2,1]
m2=dist.isetting[1,isetting]; s2=dist.isetting[2,isetting]

y1=rnorm(n[1],m1,s1)
A=cbind(rbinom(n[2],1,lambda), rnorm(n[2], m2, s2))
B=cbind(1-A[,1], rnorm(n[1], m1, s1))
y2=A[,1]*A[,2]+B[,1]*B[,2]

h1 = function(a) a; h2 = function(a) a^2
}

if (distn=="beta"){
dist.isetting=beta.settings
alpha1=dist.isetting[1,1];  beta1=dist.isetting[2,1]
alpha2=dist.isetting[1,isetting];  beta2=dist.isetting[2,isetting]

y1=rbeta(n[1], alpha1, beta1)
A=cbind(rbinom(n[2],1,lambda),rbeta(n[2], alpha2, beta2))
B=cbind(1-A[,1], rbeta(n[1], alpha1, beta1))
y2=A[,1]*A[,2]+B[,1]*B[,2]

h1=function(a) log(a); h2=function(a) log(1-a)
}

if (distn=="gamma"){
dist.isetting=gamma.settings
shape1=dist.isetting[1,1]; scale1=dist.isetting[2,1]
shape2=dist.isetting[1,isetting];scale2=dist.isetting[2,isetting]

y1=rgamma(n[1], shape=shape1, scale=scale1)
A=cbind(rbinom(n[2],1,lambda),rgamma(n[2], shape=shape2, scale=scale2))
B=cbind(1-A[,1], rgamma(n[1], shape=shape1, scale=scale1))
y2=A[,1]*A[,2]+B[,1]*B[,2]

h1=function(a) a; h2=function(a) log(a)
}



if (distn=="t"){
dist.isetting=t.settings
ncp1=dist.isetting[1,1]; df1=dist.isetting[2,1]
ncp2=dist.isetting[1,isetting]; df2=dist.isetting[2,isetting]

y1=rt(n[1], df=df1, ncp=ncp1)

A=cbind(rbinom(n[2],1,lambda), (rt(n[2], df=df2, ncp=ncp2)))
B=cbind(1-A[,1], (rt(n[1], df=df1, ncp=ncp1)))
y2=A[,1]*A[,2]+B[,1]*B[,2]

h1 = function(a) a; h2 = function(a) a^2
}

if (distn=="nbinom"){
dist.isetting=nbinom.settings
size1=dist.isetting[1,1]; p1=dist.isetting[2,1]
size2=dist.isetting[1,isetting]; p2=dist.isetting[2,isetting]

y1=rnbinom(n[1], size1, p1)+0.0001
A=cbind(rbinom(n[2],1,lambda),rnbinom(n[2], size2, p2))
B=cbind(1-A[,1], rnbinom(n[1], size1, p1))
y2=A[,1]*A[,2]+B[,1]*B[,2]+0.0001

h1=function(a) a; h2=function(a) log(a)
}

##############################################


z=c(y1,y2)
h1y1=h1(y1); h1y2=h1(y2)
h2y1=h2(y1); h2y2=h2(y2)

##############################################################
nploglik=function(mypar)
{
res <- .C("mylikC", as.integer(n[1]), as.integer(n[2]),as.double(h1y1),as.double(h2y1), as.double(h1y2), as.double(h2y2), as.double(mypar),result=double(1))
return(-res[["result"]]/nn)
}


#######modified PLRT###################################################
nmploglik=function(t)
{
penalty=cc*log(t[4])
return(nploglik(t)-penalty)
}

for(i in 2:lamda_length)
{
for(j in 1:alpha_length)
{
nmploglik_beta=function(t)
{
nmploglik(c(alpha_ini[j],t,lamda_ini[i]))
}
op.nmplrt=optim(beta.init, nmploglik_beta, hessian=TRUE, method = "Nelder-Mead", control = list(fnscale=1,maxit=1000))	
mplrt_ini[i,j]=2*(-op.nmplrt$value+nmploglik(c(0,0,0,1)))
}
}

mplrt.TS=max(mplrt_ini)



########ETMM EM  #################

TS=hat_lamda=hat_alpha=array(NA, c(lamda_length-1, alpha_length))
hat_beta=array(NA, c(2, lamda_length-1, alpha_length))


alpha.init=0

for (i in 2: lamda_length)
{
for (j in 1: alpha_length)
{

lamda.temp=lamda_ini[i]
alpha.temp=alpha_ini[j]

for (k in 1:K){
my.lamda=lamda.temp
my.alpha=alpha.temp

nmploglik_EM_beta=function(t, alpha, lamda)
{
nmploglik(c(alpha,t,lamda))
}

op=optim(beta.init, nmploglik_EM_beta, alpha=my.alpha, lamda=my.lamda,hessian=TRUE, method = "Nelder-Mead", control = list(fnscale=1,maxit=1000))	
my.beta=op$par
upper=my.lamda*exp(my.alpha+my.beta[1]*h1y2+my.beta[2]*h2y2)
lower=1-my.lamda+my.lamda*exp(my.alpha+my.beta[1]*h1y2+my.beta[2]*h2y2)

w=upper/lower

if(any(is.na(w))==FALSE){
my.lamda=(sum(w)+1)/(length(y2)+1)

op=optim(beta.init, nmploglik_EM_beta, alpha=my.alpha,lamda=my.lamda,hessian=TRUE, method = "Nelder-Mead", control = list(fnscale=1,maxit=1000))	
my.beta=op$par
alpha.temp=my.alpha
lamda.temp=my.lamda}

if(any(is.na(w))==TRUE){
my.beta=my.beta
alpha.temp=my.alpha
lamda.temp=my.lamda}

}

hat_lamda[i-1, j]=my.lamda
hat_alpha[i-1, j]=my.alpha
hat_beta[,i-1, j]=my.beta
TS[i-1, j]=2*(-nmploglik(c(my.alpha, my.beta, my.lamda))+nmploglik(c(0,0,0,1)))
}
}
mplrt_EM.TS=max(TS)


##########modified empirical LRT (Liu et al,2012)######
nloglik=function(t)
{
epis=sum((t[3]*exp(t[1]+t[2]*y2))/(1-t[3]+t[3]*exp(t[1]+t[2]*y2)))/nn
second=sum(log(1-t[3]+t[3]*exp(t[1]+t[2]*y2)))
first=sum(log(1+epis*(exp(t[1]+t[2]*z)-1)))
penalty=log(t[3])
return(-second+first-penalty)
}

for(i in 2:lamda_length)
{
nloglik_lamda=function(t)
{
nloglik(c(t[1],t[2],lamda_ini[i]))
}
ini=c(0,0)
op=optim(ini,nloglik_lamda,method = "BFGS")
lrt_ini[i]=2*(-op$value+nloglik(c(0,0,1)))
}
liu.TS=max(lrt_ini)



########score test (Qin and Liang,2011)###########
nloglik_lamda=function(t)
{
nloglik(c(t[1],t[2],1))
}
ini=c(0,0)
op=optim(ini,nloglik_lamda,method = "BFGS")
est=op$par
qin.TS=sum(exp(est[1]+est[2]*y2)-1)/(1+n[2]/n[1])


#######T test##########
t.TS=t.test(y1,y2)$statistic^2

########wilconxon test#####
wilcox.p=wilcox.test(y1,y2)$p.value

###### Logisti wald #################
logistic.TS=logistic.wald1(y1, y2, h1, h2)

return(list(mplrt_EM.TS=mplrt_EM.TS, qin.TS=qin.TS, liu.TS=liu.TS, t.TS=t.TS, wilcox.p=wilcox.p, logistic.TS=logistic.TS))
}









