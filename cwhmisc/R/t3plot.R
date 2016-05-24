T3plot <- function(x,lab=paste("T3 plot of ",deparse(substitute(x))),legend.pos="bottom", cex=0.6, ...) {
  T3 <- function(x,v) {
  # calculation of the 3rd derivative of log(m <- n(v)) for
  # nsimul simulated series, each with ndata observations
  # x: vector with observations (nsimul simulations)
    v <- cbind(v)
    n <- length(x)
    ndata <- n
    m <- nrow(v)
    xx <- matrix(rep(x,m),ncol=m,byrow=FALSE)
    vv <- matrix(rep(v,n),ncol=m,byrow=TRUE)
    sumvec <- rbind((1:ndata)/(1:ndata))
    m0 <- sumvec%*%(1/ndata*exp(xx*vv))
    m1 <- 1/ndata*sumvec%*%(xx*exp(xx*vv))
    m2 <- 1/ndata*sumvec%*%(xx^2*exp(xx*vv))
    m3 <- 1/ndata*sumvec%*%(xx^3*exp(xx*vv))
    (m3*m0-3*m2*m1+2*m1^3/m0)/m0^2
  }

  x1 <- x[is.na(x)==FALSE]
  ndata <- length(x1)
  vmax <- 1
  delta <- 0.05
  nsteps <- trunc(vmax/delta)
  sqrtn <- sqrt(ndata)
  
  #--- values in the interval [-vmax,vmax] for which 3rd derivative
  #--- will be calculated
  
  v <- -vmax+(0:(nsteps-1))*delta
  v <- c(v,0,-v[nsteps:1])
  
  #--- standardization of data
  
  x1 <- (x1-mean(x1))/sqrt(var(x1))
  
  #--- calculation of 3rd derivative Tsimul
  
  Tsimul <- sqrtn*T3(x1,v)
  
  SD <- c(9.61361441746065, 8.64515922634257, 7.79296148487048,
  7.04237226412238, 6.38080584272215,
  5.79744863963072, 5.28301418419464, 4.82953690371646,
  4.43019863437543, 4.07918264344632,
  3.77155061867714, 3.50313856389591, 3.27046787778314,
  3.07066814694318, 2.90140845110933,
  2.76083439189198, 2.64750877045021, 2.56035496977797,
  2.49860362758496, 2.46174487159861,
  2.44948974278318, 2.46174487159861, 2.49860362758496,
  2.56035496977797, 2.64750877045021,
  2.76083439189198, 2.90140845110933, 3.07066814694318,
  3.27046787778314, 3.50313856389591,
  3.77155061867714, 4.07918264344632, 4.43019863437543,
  4.82953690371646, 5.28301418419464,
  5.79744863963072, 6.38080584272215, 7.04237226412238,
  7.79296148487048, 8.64515922634257,
  9.61361441746065)
  
  
  a1 <- 0.7168311; a2 <- -2.327602; a3 <- 3.688362
  
  m <- a1+a2/sqrt(ndata)+a3/ndata
  
  # confidence limits
  
  z1 <- qnorm(0.99)
  Tz1 <- (m+z1)*SD
  
  z5 <- qnorm(0.95)
  Tz5 <- (m+z5)*SD
  
  
  # plot of T3 and confidence limits
  ym <- min(Tsimul,-Tz1)
  plot(v,Tsimul,ylim=c(ym,max(Tsimul,Tz1)),type="l",
  xlab="t",ylab="T3",main=lab,cex=0.8)
  lines(v,Tz1,lty=2,col=2)
  lines(v,-Tz1,lty=2,col=2)
  lines(v,Tz5,lty=4,col=3)
  lines(v,-Tz5,lty=4,col=3)
  legend(x=legend.pos,c("1%", "5%"), col = c(2,3), text.col= "black", lty = c(2,4), merge = TRUE, bg='white', cex=cex, ...)
}
