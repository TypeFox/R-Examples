wscv.vst <- function(alpha, theta, rho, k, n){

    z.alpha2 <- qnorm(1-alpha/2)

    rho.star <- rho / (1-rho)
    c <- 2 * (1-1/k) * (1+k*rho.star)
    f.theta <- sqrt((k-1)/2) *
      log(((1 + c*theta^2)^(0.5) - 1) / ((1 + c*theta^2)^(0.5) + 1))
    xi1 <- f.theta + z.alpha2/sqrt(n)
    xi2 <- f.theta - z.alpha2/sqrt(n)
    ci.vst.upper <- (2 * exp(xi1*(2*(k-1))^(-0.5))) /
                    (c^(0.5) * (1-exp(xi1*(2/(k-1))^(0.5))))
    ci.vst.lower <- (2 * exp(xi2*(2*(k-1))^(-0.5))) /
                    (c^(0.5) * (1-exp(xi2*(2/(k-1))^(0.5))))

    list(lbound=ci.vst.lower, ubound=ci.vst.upper)
}
#--------------------------------------------------------------------------
wscv.delta <- function(alpha, theta, rho, k, n){

    z.alpha2 <- qnorm(1-alpha/2)

    var.theta <- theta^4/(k*n) * (1 + k*rho/(1-rho)) + theta^2/(2*n*(k-1))
    sd.theta <- sqrt(var.theta)
    ci.theta <- c(theta - z.alpha2*sd.theta, theta + z.alpha2*sd.theta)

    list(lbound=ci.theta[1], ubound=ci.theta[2])
}

#---------------------------------------------------------------------------
agree.wscv <- function(ratings, conf.level=0.95, method=c("vst", "delta"), NAaction=c("fail", "omit")){

    if(!is.matrix(ratings) || ncol(ratings) < 2 || nrow(ratings) < 2)
      stop("'ratings' has to be a matrix of at least two columns and two rows.")

    na <- match.arg(NAaction)
    ratings <- switch(na,
                      fail = na.fail(ratings),
                      omit = na.omit(ratings))
    if(!is.matrix(ratings) || ncol(ratings) < 2|| nrow(ratings) < 2)
      stop("'ratings' has to be a matrix of at least two columns and two rows after removing missing values.")

    
    method <- match.arg(method)

    alpha <- 1 - conf.level
    k <- ncol(ratings)		
    n <- nrow(ratings)
    N <-	 n * k
    
    bar.x <- rowMeans(ratings)
    mu.hat <- mean(bar.x)

    sse <- sum((ratings - bar.x)^2)
    mse <- sse / (N - n)
    msb <- sum((bar.x - mu.hat)^2) * k / (n-1)

    theta <- sqrt(mse) / mu.hat

    rho <- ((n-1)*msb - n*mse) / ((n-1)*msb + n*(k-1)*mse)
    
    CI <- switch(method,
                 vst = wscv.vst(alpha, theta, rho, k, n),
                 delta = wscv.delta(alpha, theta, rho, k, n))
    list(value=theta, lbound=CI[[1]], ubound=CI[[2]])
}
