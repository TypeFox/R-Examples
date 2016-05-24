E.Quantile <- function(y, Qn, Pik) {
y<-as.data.frame(y)
Total<-rep(NA,dim(y)[2])

  if (missing(Pik))
  Pik <- rep(1, dim(y)[1])
  if (any(Pik < 0))
  stop("Probabilities must be positive.")

w <- 1/Pik
n <- length(w)

for(i in 1:dim(y)[2]){

ord <- order(y[,i])
x <- y[ord,i]
w <- w[ord]
wcum <- cumsum(w)
wsum <- wcum[n]
wper <- wsum*Qn
lows <- (wcum <= wper)
k <- sum(lows)
  if (k!=0 && k!=n){
      wlow <- wcum[k]
      whigh <- wsum - wlow
          if (whigh > wper)
          Total[i]<-x[k+1]
          else
          Total[i]<-(wlow*x[k] + whigh*x[k+1]) / wsum
  }
  if (k == 0) {
      Total[i] <- x[1]
  }
  if (k == n) {
      Total[i] <- x[n]
  }
}
return(Total)
}
