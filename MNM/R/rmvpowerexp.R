`rmvpowerexp` <-
function(n, Location = rep(0, nrow(Scatter)), Scatter = diag(length(Location)), Beta=1) 
        {
        p <- length(Location)
        if (!isSymmetric(Scatter, tol = sqrt(.Machine$double.eps))) {stop("Scatter must be a symmetric matrix")}
        if (p != nrow(Scatter)) {stop("Location and Scatter have non-conforming size")}

        ev <- eigen(Scatter, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {warning("Scatter is numerically not positive definite")}
        ScatterSqrt <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
        
        radius <- (rgamma(n, shape = p/(2*Beta), scale = 1/2))^(1/(2*Beta))
        
        un <- runifsphere(n=n,p=p)
        mvpowerexp <- radius * un %*% ScatterSqrt
        mvpowerexp <- sweep(mvpowerexp, 2, Location, "+")

        return(mvpowerexp)
        }
