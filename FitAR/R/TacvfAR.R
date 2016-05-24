`TacvfAR` <-
function(phi, lag.max = 20)
{
    p<-length(phi)
    maxlagp1<-lag.max+1
    if(p == 0)
        if(lag.max >= 0) 
            return(c(1, numeric(lag.max)))
        else 
            stop("maxlag invalid")
    r <- p + 1
    b <- numeric(r)
    C <- 1
    phiStar <- numeric(3 * r)
    phiStar[r] <- -1
    phiStar[r + 1:p]<-phi
    a <- matrix(numeric(r^2), ncol = r)
    b[1] <- 1
    for(i in 1:r)
        for(j in 1:r)
           if(j == 1)
                  a[i, j] <- phiStar[r + i - 1]  
           else 
                  a[i, j] <- phiStar[r + i - j] + phiStar[r + i + j - 2]
    g <- solve(a,  -b)
    if(length(g) < maxlagp1) {
        g <- c(g, numeric(maxlagp1 - r))
        for(i in (r + 1):maxlagp1) 
            g[i] <- phi %*% g[i - 1:p]
    }
    g[1:maxlagp1]
}

