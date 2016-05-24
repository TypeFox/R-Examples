cpoint <- function(x, mu, sigma){
 DELTA <- deltat(x)
 n <- length(x)
 Z <- NULL
 if(!missing(mu) && !missing(sigma)){    
  Z <- (diff(x) - mu(x[1:(n-1)])*DELTA)/(sqrt(DELTA)*sigma(x[1:(n-1)]))
 } else {
	bw <- n^(-1/5) * sd(x)
    y <- sapply(x[1:(n-1)], function(xval) {
        tmp <- dnorm(xval, x[1:(n - 1)], bw)
        sum(tmp * diff(x))/(DELTA * sum(tmp))
    })
   Z <- diff(x)/sqrt(DELTA) - y*sqrt(DELTA)    
 }
 
 lenZ <- length(Z)	
 Sn <- cumsum(Z^2)
 S <- sum(Z^2)
 D <- abs(1:lenZ/lenZ - Sn/S)
 k0 <- which(D==max(D))[1]
 return(list(k0=k0+1, tau0=time(x)[k0+1], 
   theta1=sqrt(Sn[k0]/k0), theta2=sqrt((S-Sn[k0])/(lenZ-k0))))
}
