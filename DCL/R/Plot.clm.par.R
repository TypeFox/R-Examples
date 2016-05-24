Plot.clm.par<-function(clm.par)
{ 
  # The alpha's
  alpha<-clm.par$alpha
  m<-length(alpha)
  # The beta's
  beta<-clm.par$beta
  # The development factors
  Fj<-clm.par$Fj
  par(mfrow=c(3,1))
  par(mar=c(3,1.5,1.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(1.5,0.5,0))
  # c(bottom, left, top, right)
  par(cex.main=1)
  plot(1:m,alpha,type='b' ,pch=19,col=1,cex=0.8,xlab='underwriting period',
       ylab='', main='CL underwriting parameters',xaxt='n')
  axis(1:m,as.character(1:m))
  plot(1:m,beta,type='b' ,pch=19,lty=1,col=1,cex=0.8,xlab='development period',
       ylab='', main='CL development parameters',xaxt='n')
  axis(1:m,as.character(0:(m-1)))
  plot(1:(m-1),Fj,type='b' ,pch=19,lty=1,col=1,cex=0.8,xlab='development period',
       ylab='', main='Development factors',xaxt='n')
  axis(1:m,as.character(1:(m-1)))
  
  par(mfrow=c(1,1)) 
}  
