.householder <-
function(x){
        m <- length(x)
        alpha <- sqrt(drop(crossprod(x)))
        e <- c(1,rep(0,m-1))
        u <- x - alpha*e
        v <- u/sqrt(drop(crossprod(u)))
        diag(m) - 2*v%*%t(v)
    }
