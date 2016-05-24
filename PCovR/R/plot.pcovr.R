plot.pcovr <-
  function(x,cpal=NULL,lpal=NULL, ...){
    modsel <- x$modsel
    alpha <- x$alpha
    R <- x$R
    vec <- x$Rvalues
    a <- x$Alphavalues
    if (is.null(cpal)){cpal <- 1:length(vec)}
    if (is.null(lpal)){lpal <- 1:length(vec)}
    
    if (modsel=="sim"){
      Qy2 <- x$Qy2
      plot(a,Qy2[,1],xlab="Weighting Parameter",ylab="Qy2",col=cpal[1],lty=lpal[1],type="l",ylim=c(0,1))
      for (r in 1:length(vec)){
        points(a,Qy2[,r],col=cpal[r],type="l",lty=lpal[r])
      }
      points(alpha,Qy2[which(a==alpha),which(vec==R)],pch=2,col=cpal[vec==R])
      text(alpha,Qy2[which(a==alpha),which(vec==R)]+.05,"optimal",col=cpal[vec==R])
      legend("topleft",as.character(vec),title="Number of Components",col=cpal,lty=lpal)
      
    } else if (modsel=="seqRcv"){
      Qy2 <- x$Qy2
      plot(vec,Qy2,xlab="Number of Components",ylab="Qy2",col=cpal[1],lty=lpal[1],type="l",ylim=c(0,1))
      points(R,Qy2[vec==R],pch=2,col=cpal[1])
      text(R,Qy2[vec==R]-.05,"optimal")
      legend("topleft",c("Qy2"),col=cpal[1],lty=lpal[1])
      
    } else if (modsel=="seq"){
      VAF <- x$VAFsum
      plot(vec,VAF,xlab="Number of Components",ylab="VAFsum",col=cpal[1],lty=lpal[1],type="l",ylim=c(0,1))
      points(R,VAF[vec==R],pch=2,col=cpal[1])
      text(R,VAF[vec==R]-.05,"optimal")
      legend("topleft",c("VAFsum"),col=cpal[1],lty=lpal[1])
      
    } else if (modsel=="seqAcv"){
      par(mfrow = c(2,1))
      VAF <- x$VAFsum
      Qy2 <- data.matrix(x$Qy2)
      plot(vec,VAF,xlab="Number of Components",ylab="VAFsum",col=cpal[1],lty=lpal[1],type="l",ylim=c(0,1))
      points(R,VAF[vec==R],pch=2,col=cpal[1])
      text(R,VAF[vec==R]-.05,"optimal")
      legend("topleft",c("VAFsum"),col=cpal[1],lty=lpal[1])
      plot(a,Qy2,xlab="Weighting parameter",ylab="Qy2",col=cpal[1],lty=lpal[1],type="l",ylim=c(0,1))
      points(alpha,t(Qy2)[a==alpha],pch=2,col=cpal[1])
      text(alpha,Qy2[a==alpha]-.05,"optimal")
      legend("topleft",c("Qy2"),col=cpal[1],lty=lpal[1])
    } 
  }