siegel.tukey <-
function(x,y)  {
  m <- length(x);n <- length(y)
  N <- m+n
  an <- function(N){
    TEMP <- NULL
    for(i in 1:N){
      if(i<=N/2){
        if(i%%2==0) TEMP[i] <- 2*i else TEMP[i] <- 2*i-1
      }
      if(i>N/2){
        if(i%%2==0) TEMP[i] <- 2*(N-i)+2 else TEMP[i] <- 2*(N-i)+1
      }
    }
    return(TEMP)
  }
  z <- function(x,y){
    TEMP2 <- cbind(c(x,y),c(rep(1,length(x)),rep(0,length(y))))
    return(TEMP2[order(TEMP2[,1]),2])
  }
  sn <- sum(z(x,y)*an(N))
  zz <- (sn-m*(N+1)/2)/sqrt(m*n*(N+1)/12)
  return(c(zz,pnorm(zz)))
}
