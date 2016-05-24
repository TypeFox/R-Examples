lssd <- function(data,k1,p = 2,q = 0,interval = c(0.001,0.999)) {
    
    n <- length(data)
    k <- 1:k1
    nk <- n/k
    kp1 <- 1/(k^p)
    logk <- log(k)
    
    sdrunning <- function(kscale) {
        sd(running(data,fun = sum,width = kscale,by = kscale))
    }
    
    logs <- log(sapply(k,sdrunning))
    a1 <- sum(kp1)

    logck <- function(H) {0.5*log((nk - nk^(2*H - 1))/(nk - 0.5))}
    
    g1 <- function(H) {H * logk + logck(H) - logs}
    g2 <- function(H) {- (sum(kp1 * g1(H)))^2 + a1 * sum(kp1 * (g1(H))^2)}
    g3 <- function(H) {g2(H) + (H^(q+1))/(q+1)}
    
    Hest <- ifelse(q == 0,optimize(g2,interval = interval),
    optimize(g3,interval = interval))[[1]]
    
    sigmaest <- exp((sum(logs*kp1) - Hest * sum(logk*kp1) -
    sum(logck(Hest)*kp1))/a1)
    
    return(setNames(c(sigmaest,Hest),c("sigma_estimate","H_estimate")))
}