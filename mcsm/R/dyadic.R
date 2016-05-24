dyadic=function(N=10^4,q=4){
# A dyadic improvement for a toy example

  h=function(x){(cos(50*x)+sin(20*x))^2}
  uref=runif(N)
  x=h(uref)
  estx=cumsum(x)/(1:N)

  resid=uref%%2^(-q)
  simx=matrix(resid,ncol=2^q,nrow=N)
  simx[,2^(q-1)+1:2^1]=2^(-q)-simx[,2^(q-1)+1:2^1]
  for (i in 1:2^q) simx[,i]=simx[,i]+(i-1)*2^(-q)
  xsym=h(simx)
  estint=cumsum(apply(xsym,1,mean))/(1:N)

  q=2*q
  resid=uref%%2^(-q)
  simx=matrix(resid,ncol=2^q,nrow=N)
  simx[,2^(q-1)+1:2^1]=2^(-q)-simx[,2^(q-1)+1:2^1]
  for (i in 1:2^q) simx[,i]=simx[,i]+(i-1)*2^(-q)
  xsym=h(simx)
  estunt=cumsum(apply(xsym,1,mean))/(1:10^4)

  plot(estx,ty="l",xlab="Iterations", ylab="",lwd=2,ylim=c(.8,1.))
  lines(estint,lwd=2,lty=2,col="grey30")
  lines(estunt,lwd=2,col="sienna")
  }
