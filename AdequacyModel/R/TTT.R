TTT <-
function(x,lwd=2,lty=2,col="black",grid=TRUE,...)
{
  
  n <- length(x)
  y <- order(x)
  sortx <- x[y]
  Trn <- rep(NA, times = n)
  r   <- rep(NA, times = n)
  Trn[1] <- n*sortx[1]
  r[1] <- 1/n
  
  for(i in 2:n)
  {
    Trn[i] <- Trn[i-1]+(n-i+1)*(sortx[i]-sortx[i-1])
    r[i]   <- i/n
  }
  
  plot(r,Trn/Trn[n],xlab="i/n",ylab="T(i/n)",xlim=c(0,1), ylim=c(0,1), main="",type="l",lwd=lwd, lty=1,col=col,...)
  lines(c(0,1),c(0,1),lty=lty,lwd=1,...)
  if(grid=="TRUE") grid()
}
