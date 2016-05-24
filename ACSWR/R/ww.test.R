ww.test <-
function(x,y)  {
  runfunction <- function(x,y){
    xind <- rep(1,length(x))
    yind <- rep(2,length(y))
    xy <- c(x,y); xyind <- c(xind,yind);grand <- cbind(xy,xyind)
    grand <- grand[rank(grand[,1]),]
    num_of_runs <- sum(diff(grand[,2])!=0)+1
    # return(num_of_runs)
  }
  m <- length(x); n <- length(y)
  mu0 <- 1+2*m*n/(m+n)
  var0 <- 2*m*n*(2*m*n-m-m)/((m+n-1)*(m+n)^2)
  test.statistic <- (runfunction(x,y)-mu0)/sqrt(var0)
  return(2*(1-pnorm(abs(test.statistic))))
}
