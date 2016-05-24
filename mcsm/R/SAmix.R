SAmix=function(x,tolerance=10^(-4),factor=1){
#SA version

like=function(mu){
  -sum(log((.25*dnorm(da-mu[1])+.75*dnorm(da-mu[2]))))
  }

temp=scale=iter=dif=1
the=matrix(x,ncol=2)
curlike=hval=like(x)

while (dif>tolerance){

  prop=the[iter,]+rnorm(2)*scale[iter]
  if ((max(-prop)>2)||(max(prop)>5)||(temp[iter]*log(runif(1))>-like(prop)+curlike))
      prop=the[iter,]

  curlike=like(prop);hval=c(hval,curlike);the=rbind(the,prop)
  iter=iter+1;temp=c(temp,1/log(iter+1))
  ace=length(unique(the[(iter/2):iter,1]))
  if (ace==1) factor=factor/10
  if (2*ace>iter) factor=factor*10
  scale=c(scale,max(2,factor*sqrt(temp[iter])))
  dif=(iter<100)+(ace<2)+(max(hval)-max(hval[1:(iter/2)]))
  }

  list(theta=the,like=hval,ite=iter)
}
