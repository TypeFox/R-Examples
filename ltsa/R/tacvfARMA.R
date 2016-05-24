`tacvfARMA` <-
function(phi = numeric(0), theta = numeric(0), maxLag = 1, sigma2=1)
{
stopifnot(maxLag >= 0, sigma2>0)
#check for stationarity
    p <- length(phi)
    if(p>0){
        pi=numeric(p)
        phik <- phi
            for (k in 1:p){
                L <- p+1-k
                a <- phik[L]
                pi[L+1-k] <- a
                phikp1 <- phik[-L]
                if(is.na(a) || abs(a)==1) 
                    stop("error: roots outside unit circle -- nonstationary/noncausal model")
            phik <- (phikp1+a*rev(phikp1))/(1-a^2)
            }
    if (!all(abs(pi)<1)) 
        stop("error: roots outside unit circle -- nonstationary/noncausal model")
    }
#model is stationary, compute acvf
    q <- length(theta)
    if(max(p, q) == 0) 
        return(c(sigma2, numeric(maxLagp1))) #white noise case
    maxLagp1 <- maxLag + 1
    r <- max(p, q) + 1
    b <- numeric(r)
    C <- numeric(q + 1)
    C[1] <- 1
    theta2 <- c(-1, theta)
    phi2 <- numeric(3 * r)
    phi2[r] <- -1
    if(p > 0) 
        phi2[r + 1:p] <- phi
    if(q > 0) 
        for(k in 1:q) {
            C[k + 1] <-  - theta[k]
            if(p > 0) 
                for(i in 1:min(p, k)) 
                  C[k + 1] <- C[k + 1] + phi[i] * C[k + 1 - i]
        }
    for(k in 0:q) 
        for(i in k:q) 
            b[k + 1] <- b[k + 1] - theta2[i + 1] * C[i - k + 1]
    if(p == 0) 
        g <- c(b, numeric(maxLagp1))[1:maxLagp1]
    else {
        a <- matrix(numeric(r^2), ncol = r)
        for(i in 1:r) 
            for(j in 1:r) {
                if(j == 1) 
                  a[i, j] <- phi2[r + i - 1]
                else 
                  a[i, j] <- phi2[r + i - j] + phi2[r + i + j - 2]
            }
        g <- solve(a,  - b)
        if(length(g) <= maxLag) {
            g <- c(g, numeric(maxLag - r))
            for(i in (r + 1):maxLagp1) 
                g[i] <- phi %*% g[i - 1:p]
          } 
    }
    sigma2*g[1:maxLagp1]
}

