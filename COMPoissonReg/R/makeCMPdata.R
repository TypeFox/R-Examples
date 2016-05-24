makeCMPdata <- function(x,beta,nu){

# Given x, create full x matrix
  #create vector of ones
     if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x))

  #create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
     newx <- cbind(onevec,x)
     xmat <- as.matrix(newx)

# Create lambda vector
  lambda <- exp(xmat %*% beta)


# generate uniform distribution values
  unifvals <- runif(n=length(lambda))

# Create space for simulated y's
  keepy <- rep(0,length(lambda))


 for (i in 1:length(unifvals)){
     # start counter for y
       y <- 0
     # Compute Z-inverse.  This equals P(Y=0).
       zinv <- 1/computez.lambdaest(lambda[i],nu,max=200)
       py <- zinv
       while (py < unifvals[i]){
           y <- y+1
           py <- py + ((lambda[i]^y)/((factorial(y))^nu))*zinv
          }
     # Keep last value of y, ie where py > unifvals[i]
     keepy[i] <- y
  }

return(keepy)
}

