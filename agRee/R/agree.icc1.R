icc1.sf <- function(alpha, msb, mse, n, k){

    F0 <- msb / mse
    FU <- F0 * qf(1-0.5*alpha, n*(k-1), n-1)
    FL <- F0 / qf(1-0.5*alpha, n-1, n*(k-1))

    list(lbound=(FL-1) / (FL+k-1), ubound=(FU-1) / (FU+k-1))
}
#--------------------------------------------------------------------------
agree.icc1 <- function(ratings, conf.level=0.95, method=c("sf"), NAaction=c("fail", "omit")){

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

    rho <- (msb - mse) / (msb + (k-1)*mse)
    
    CI <- switch(method,
                 sf = icc1.sf(alpha, msb, mse, n, k))
    
    list(value=rho, lbound=CI[[1]], ubound=CI[[2]])
}
