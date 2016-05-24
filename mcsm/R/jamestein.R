jamestein=function(N=10^3,p=5){
#Plot of Monte Carlo approximation to James-Stein risks

  nor=matrix(rnorm(N*p),nrow=p)
  risk=matrix(0,ncol=100,nrow=10)
  a=seq(1,2*(p-2),le=10)
  the=sqrt(seq(0,6*p,le=100)/p)
  for (j in 1:100){
      para=rep(the[j],p)
      nornor=apply((nor+para)^2,2,sum)
      shrink=(1-outer(a,nornor,"/"))
      shrink=shrink*(shrink>0)
      for (i in 1:10){
      for (t in 1:N)
        risk[i,j]=risk[i,j]+sum((para-shrink[i,t]*(nor[,t]+para))^2)
        }}
  risk=risk/N
  par(mar=c(4,2,2,1))
  plot(sqrt(p)*the,risk[1,],ty="l",yl=c(0,p),lwd=2,xlab=expression(theta),ylab="")
  for (i in 2:10) lines(sqrt(p)*the,risk[i,],lwd=2)
}
