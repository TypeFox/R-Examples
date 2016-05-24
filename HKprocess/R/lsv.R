lsv <- function(data,k1,p = 6,q = 0,interval = c(0.001,0.999)) {
    
    n <- length(data)
    k <- 1:k1
    sdrunning <- function(kscale) {
        sd(running(data,fun = sum,width = kscale,by = kscale))
    }
    
    s <- sapply(k,sdrunning)
    d1 <- sum((s^4)/(k^p))
    
    # eq.22 in Tyralis and Koutsoyiannis (2011)
    g2 <- function(H) {
        ckH <- ((n/k) - (n/k)^(2*H - 1))/(n/k - 1)
        a1H <- sum(((ckH^2) * (k^(4 * H)))/(k^p))
        a2H <- sum((ckH * (k^(2 * H)) * (s^2))/(k^p))
        d1 - (a2H^2)/a1H
    }
    
    # eq.24 in Tyralis and Koutsoyiannis (2011)
    g3 <- function(H) {g2(H) + (H^(q+1))/(q+1)}
    
    Hest <- ifelse(q == 0,optimize(g2,interval = interval),
    optimize(g3,interval = interval))[[1]]
    
    ckHest <- ((n/k) - (n/k)^(2*Hest - 1))/(n/k - 1)
    a1Hest <- sum(((ckHest^2) * (k^(4 * Hest)))/(k^p))
    a2Hest <- sum((ckHest * (k^(2 * Hest)) * (s^2))/(k^p))
    sigmaest <- sqrt(a2Hest/a1Hest)
    return(setNames(c(sigmaest,Hest),c("sigma_estimate","H_estimate")))
}