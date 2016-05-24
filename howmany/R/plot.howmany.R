"plot.howmany" <-
function(x,...)
  {
    m <- x$m
    m0 <- length(x$pvalues)
    par(mfrow=c(2,2))
    
    plot(x$pvalues,1:m0,main="distribution of  p-values",type="l",xlab="p",ylab="number of rejected hypotheses")
    lines(x$pvalues,x$boundingfunction,lty=3)
    lines(c(0,1),c(0,m),lty=2)

    plot(x$pvalues,(1:m0)-x$pvalues*m,type="l",main="excess p-values",xlab="p",ylab="excess p-values")
    abline(h=0,lty=2)
    lines(x$pvalues,x$boundingfunction-x$pvalues*m,lty=3)

    plot(1:m0,x$lowerbound,xlab="number of rejected hypotheses",ylab="number of correct rejections",type="l",main="correct rejections\n (simultaneous lower bound)")

    plot(1:m0,x$lowerbound/(1:m0),xlab="number of rejected hypotheses",ylab="proportion of correct rejections",type="l",main="proportion of correct rejections\n (simultaneous lower bound)")

  }

