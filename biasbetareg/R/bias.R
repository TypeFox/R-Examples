bias=function(fit){

mu= predict(fit,type="response")
phi=predict(fit,type="precision")
n=length(fit$y)
k=length(fit$coefficients$mean)
h=length(fit$coefficients$precision)

a=trigamma((1-mu)*phi)+trigamma((mu)*(phi))
b=(trigamma((1-mu)*phi))*((1-mu)^2)+(trigamma(mu*phi))*((mu)^2)-trigamma(phi)
c=psigamma((1-mu)*phi,deriv=2)+psigamma(mu*phi,deriv=2)
d=(psigamma((1-mu)*phi,deriv=2))*((1-mu)^2)-(psigamma(mu*phi,deriv=2))*((mu)^2)
e=(psigamma(phi,deriv=2))-(psigamma(mu*phi,deriv=2))*((mu)^3)-(psigamma((1-mu)*phi,deriv=2))*((1-mu)^3)


if(fit$link$mean$name=="logit"){

du=(mu*(1-mu))
du2=(mu*(1-mu)*(1-2*mu))

}
else if(fit$link$mean$name=="probit") {

du=(1/(sqrt(2*pi)*exp(-1/2*((pnorm(mu)^-1)^2))))
du2=(-(pnorm(mu)^-1)/(sqrt(2*pi)*exp(-1/2*((pnorm(mu)^-1)^2))))

}
else if(fit$link$mean$name=="cloglog"){

du=(-log(1-mu)*(1-mu))
du2=(-(1-mu)*log(1-mu)*(1+log(1-mu)))

}



if(fit$link$precision$name=="identity"){

df=1
df2=0

}
else if(fit$link$precision$name=="log") {

df=phi
df2=phi

}
else if(fit$link$precision$name=="sqrt"){

df=(2*sqrt(phi))
df2=2

}


Wtt=diag((b*(df^2)),n,n)
Wbb=diag((((phi)^2)*a*(du^2)),n,n)
Wbt=diag((phi*((mu*a)-trigamma((1-mu)*phi))*du*df),n,n)

M1=diag((((phi^2)/2)*(phi*c*(du^3)-a*du*du2)),n,n)
M2=diag((((phi^2)/2)*(mu*c-psigamma((1-mu)*phi,deriv=2))*(du^2)*df+phi*(trigamma((1-mu)*phi)-a*mu)*(du2)*df),n,n)
M3=diag(((-phi/2)*((2*a+phi*psigamma((1-mu)*phi,deriv=2)-phi*mu*c)*(du^2)*df+(trigamma((1-mu)*phi)-mu*a)*(du2)*df)),n,n)
M4=diag(((1/2)*((d*phi+2*trigamma((1-mu)*phi)-2*mu*a)*du*(df^2)-phi*(trigamma((1-mu)*phi)-mu*a)*du*df2)),n,n)
M5=diag(((phi/2)*(d*du*(df^2)+(trigamma((1-mu)*phi)-mu*a)*du*df2)),n,n)
M6=diag(((1/2)*(e*(df^3)-b*df*df2)),n,n)

X=matrix(model.matrix(fit),n,k)
Z=matrix(model.matrix(fit,"precision"),n,h)
O1=matrix(c(rep(0,n*h)),nrow=n)
O2=matrix(c(rep(0,n*k)),nrow=n)

A=cbind(X,O1)
B=cbind(O2,Z)
P=rbind(A,B)

C=cbind(Wbb,Wbt)
D=cbind(Wbt,Wtt)
W=rbind(C,D)

p=t(P)
K=p%*%W%*%P
K1=solve(K) 

Kbb=K1[1:k,1:k]
Kbt=K1[1:k,(k+1):(k+h)]
Ktt=K1[(k+1):(k+h),(k+1):(k+h)]

Xtil=X
Ztil=Z

pbb=(Xtil)%*%(Kbb)%*%(t(Xtil))
pbt=(Xtil)%*%(Kbt)%*%(t(Ztil))
ptt=(Ztil)%*%(Ktt)%*%(t(Ztil))

Pbb=matrix(diag(pbb),n,1)
Pbt=matrix(diag(pbt),n,1)
Ptt=matrix(diag(ptt),n,1)

w1=matrix (c(((M1%*%Pbb)+((M2+M3)%*%Pbt)+(M5%*%Ptt)),((M2%*%Pbb)+((M4+M5)%*%Pbt)+(M6%*%Ptt))),nrow=2*n)
W1=solve(W)
E1=W1%*%w1

B=K1%*%p%*%W%*%E1
m1=as.matrix(fit$coefficients$mean)
m2=as.matrix(fit$coefficients$precision)
coef=rbind(m1,m2)

cor=coef-B
cor
cbind(matrix(names(fit$bias),length(names(fit$bias)),1),cor)
 
}
