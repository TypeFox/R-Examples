brownian.motion.variance <-
function(n.locs, time.lag, location.error, x, y, max.lag){

  #Creating NULL vectors to store data
  T.jump <- alpha <- ztz <- ob <- loc.error.1 <- loc.error.2 <- NULL
  
  i <- 2
  while(i < n.locs){
	  if((time.lag[i]+time.lag[i-1])/2 >= max.lag) {
		  i = i + 1
	  } else {
		  ob <- c(ob, i)
		  t <- time.lag[i]+time.lag[i-1]
		  T.jump <- c(T.jump, t)
		  a <- time.lag[i] / t
		  alpha <- c(alpha, a)
		  u <- c(x[i-1], y[i-1]) + a*(c(x[i+1], y[i+1]) - c(x[i-1], y[i-1]))
		  ztz <- c(ztz, (c(x[i], y[i]) - u)%*%(c(x[i], y[i]) - u))
		  loc.error.1 <- c(loc.error.1, location.error[i-1])
		  loc.error.2 <- c(loc.error.2, location.error[i+1])
		  i <- i + 2
	  }
  }
  
  #Likelihood function for Brownian Motion variance estimation
  likelihood <- function(var){    
      v <- T.jump*alpha*(1-alpha)*var + ((1-alpha)^2)*(loc.error.1^2) + 
              (alpha^2)*(loc.error.2^2)
      l <- (1/(2*pi*v))*exp(-ztz/(2*v))
      -sum(log(l), na.rm=TRUE)
  }
  
  BMvar <- optimize(likelihood, lower=1, upper=1000000)$minimum

  return(BMvar)

}
