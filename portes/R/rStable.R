"rStable" <- 
function (n, ALPHA, BETA, GAMMA = NULL, DELTA = NULL) 
{
    if (any(ALPHA > 2)) 
        stop("Error: ALPHA is greater than 2")
    if (any(ALPHA <= 0)) 
        stop("Error: ALPHA is less than or equal to 0")
    if (any(BETA < -1)) 
        stop("Error: BETA is less than -1")
    if (any(BETA > 1)) 
        stop("Error: BETA is greater than 1")
    k <- length(ALPHA)
    if (length(BETA) !=k) 
        stop("Error: ALPHA and BETA should have the same dimension")
    theta <- matrix( pi * (runif(k*n) - 1/2),ncol=k)
    z <- matrix(-log(runif(k*n)),ncol=k)
      result <- ans <- matrix(numeric(k*n),ncol=k)
        c <- (1 + (BETA * tan(pi * ALPHA/2))^2)^(1/(2 * ALPHA))
        theta0 = (1/ALPHA) * atan(BETA * tan(pi * ALPHA/2))
     for (i in 1:k){
       if (ALPHA[i] == 1 & BETA[i] == 0) 
        result[,i] <- matrix(rcauchy(k*n),ncol=k)
       else
        result[,i] <- (c[i] * sin(ALPHA[i] * (theta[,i] + theta0[i]))/(cos(theta[,i]))^(1/ALPHA[i])) * (cos(theta[,i] - ALPHA[i] * (theta[,i] + theta0[i]))/z[,i])^((1 - ALPHA[i])/ALPHA[i])
       result[,i] <- result[,i] - BETA[i] * tan(ALPHA[i] * pi/2)
       if (is.null(GAMMA)) GAMMA <- rep(1,k)
       if (is.null(DELTA)) DELTA <- rep(0,k)
      ans[,i] <- result[,i] * GAMMA[i] + DELTA[i]
     }
  return(ans)
}
