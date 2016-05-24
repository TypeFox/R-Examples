EMcenso=function(repp=10){
#data generation
n=30;m=20;ybar=0.32;a=1

#True loglikelihood
like=function(the){
  dnorm(the,mean=ybar,sd=1/sqrt(m),log=T)+(n-m)*pnorm(the-a,log=T)}

#EM sequence
EMI=function(){
theta=rnorm(1,mean=ybar)
iter=nonstop=1

while (nonstop>0){
  theta=c(theta,m*ybar/n+(n-m)*(theta[iter]+
  dnorm(a-theta[iter])/pnorm(a-theta[iter]))/n)
  nonstop=(abs(diff(theta[iter:(iter+1)]))>10^(-4))
  iter=iter+1
  }
  theta
  }

#Graphs
theta=EMI()
plot(theta,like(theta),type="l",xlab=expression(theta[0]),
ylab=expression(l(theta[0])),lty=2,xlim=c(-0.5,1.8),ylim=c(-40,-10))

for (t in 1:repp){
  theta=EMI();lines(theta,like(theta),lty=2)}

curve(like,add=T,lwd=4,col="grey70")
}
