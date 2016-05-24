Plot.dcl.par<-function(dcl.par,type.inflat='DCL')
{
  par(mfrow=c(2,2))
  par(mar=c(3,1.5,1.5,1.5),oma=c(2,0.5,0.5,0.2),mgp=c(1.5,0.5,0))
  # c(bottom, left, top, right)
  par(cex.main=1)
  alpha.X<-dcl.par$alpha.X
  m<-length(alpha.X)
  beta.X<-dcl.par$beta.X
  inflat<-dcl.par$inflat
  p.delay<-dcl.par$pi.delay
  pj<-dcl.par$pj
  
  if (type.inflat=="IDCL")
  {
      inflat.DCL<-dcl.par$inflat.DCL
      alpha.I<-dcl.par$alpha.I
      yy<-range(c(alpha.I,alpha.X),na.rm=TRUE)
      plot(1:m,alpha.X,type='b' ,pch=15,lty=2,col=2,cex.axis=0.8,xlab='underwriting period',
           ylab='', ylim=yy, main='CL underwriting parameters',xaxt='n')
      axis(1,1:m,as.character(1:m))
      
      points(1:m,alpha.I,pch=15,col=4)
      lines(1:m,alpha.I,lty=1,col=4)
      legend('topleft',cex=0.8,bty='n',legend=c('incurred','paid'),col=c(4,2),
             pch=c(19,15))#,lty=c(1,2))
      grid(ny=5,nx=NA)
      plot(1:m,beta.X,type='b' ,pch=19,lty=1,col=1,cex.axis=0.8,xaxt='n',xlab='development period',
           ylab='', main='CL development parameters')
      axis(1,1:m,as.character(0:(m-1)))
      grid(ny=5,nx=NA)
      
      yy<-range(c(inflat,inflat.DCL),na.rm=TRUE)
      plot(1:m,inflat,col=4,pch=19,ylim=yy,ylab='',xlab='underwriting period',
           main='Severity inflation',cex.axis=0.8,xaxt='n')
      axis(1,1:m,as.character(1:m))
      
      grid(ny=5,nx=NA)
      lines(inflat,col=4,lty=1,lwd=1)
      
      points(inflat.DCL,col=2,pch=15)
      lines(inflat.DCL,col=2,lty=2,lwd=1)
      legend('topleft',c('IDCL','DCL'),
             col=c(4,2),pch=c(19,15),cex=0.8,bty='n')#,lty=c(1,2))
      
    } 
  
  if (type.inflat=="BDCL")
  {
    inflat.DCL<-dcl.par$inflat.DCL
    alpha.I<-dcl.par$alpha.I
        yy<-range(c(alpha.I,alpha.X),na.rm=TRUE)
        plot(1:m,alpha.X,type='b' ,pch=15,lty=2,col=2,cex.axis=0.8,xlab='underwriting period',
             ylab='', ylim=yy, main='CL underwriting parameters',xaxt='n')
        axis(1,1:m,as.character(1:m))
        
        points(1:m,alpha.I,pch=19,col=4)
        lines(1:m,alpha.I,lty=1,col=4)
        legend('topleft',cex=0.8,bty='n',legend=c('incurred','paid'),col=c(4,2),
               pch=c(19,15))#,lty=c(1,2))
        grid(ny=5,nx=NA)
        plot(1:m,beta.X,type='b' ,pch=19,lty=1,col=1,cex.axis=0.8,xaxt='n',xlab='development period',
             ylab='', main='CL development parameters')
        axis(1,1:m,as.character(0:(m-1)))
        grid(ny=5,nx=NA)
        
        yy<-range(c(inflat,inflat.DCL),na.rm=TRUE)
        plot(1:m,inflat,col=4,pch=19,ylim=yy,ylab='',xlab='underwriting period',
             main='Severity inflation',cex.axis=0.8,xaxt='n')
        axis(1,1:m,as.character(1:m))
        
        grid(ny=5,nx=NA)
        lines(inflat,col=4,lty=1,lwd=1)
        
        points(inflat.DCL,col=2,pch=15)
        lines(inflat.DCL,col=2,lty=2,lwd=1)
        legend('topleft',c('BDCL','DCL'),
               col=c(4,2),pch=c(19,15),cex=0.8,bty='n')#,lty=c(1,2))
    
  }
  
  if (type.inflat=="DCL")
  {
      plot(1:m,alpha.X,type='b' ,pch=19,col=1,cex.axis=0.8,xlab='underwriting period',
           ylab='', main='CL underwriting parameters',xaxt='n')
      axis(1,1:m,as.character(1:m))
      
      grid(ny=5,nx=NA)
      plot(1:m,beta.X,type='b' ,pch=19,lty=1,col=1,cex.axis=0.8,xaxt='n',xlab='development period',
           ylab='', main='CL development parameters')
      axis(1,1:m,as.character(0:(m-1)))
      grid(ny=5,nx=NA)
      plot(1:m,inflat,col=4,pch=19,ylab='',xlab='underwriting period',
           main='Severity inflation',cex.axis=0.8,xaxt='n')
      axis(1,1:m,as.character(1:m))
      
      grid(ny=5,nx=NA)
      lines(inflat,col=4,lty=2,lwd=1)
    
  }
  
  plot(1:m,p.delay,type='l',lwd=2,lty=1,col=4,xaxt='n',xlab='settlement delay',
       ylab='',main='Delay parameters',cex.axis=0.8)
  axis(1,1:m,as.character(0:(m-1)))
  grid(ny=5,nx=NA)
  points(1:m,p.delay,pch=18,col=4,cex=0.6)
  lines(1:m,pj,type='l',lwd=2,lty=2,col=3)
  points(1:m,pj,pch=15,col=3,cex=0.6)
  grid(ny=5,nx=NA)
  legend('topright',legend=c('general','adjusted'),col=c('blue','green'),
           bty='n',pch=c(18,15),cex=0.9)#,lty=c(1,2))
  par(mfrow=c(1,1))
}