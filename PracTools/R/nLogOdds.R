
nLogOdds <- function(moe.sw, e, alpha=0.05, pU, N=Inf){
    if (!(moe.sw==1 | moe.sw==2))
        stop("moe.sw must equal 1 or 2.\n")
    if (alpha <= 0 | alpha >= 1)
        stop("alpha must be in (0,1).\n")            
    if (sum(sapply(list(e, N, pU), is.null) != 0))
        stop("e, N, and pU cannot be NULL.\n")
    if (any(pU <= 0) | any(pU >= 1)) stop("pU must be in (0,1).\n")
    if (N <= 0) stop("N must be positve.\n")
    
    if (N == Inf) {a <- 1}
        else {a <- N/(N-1)}
            
    z <- qnorm(1 - alpha/2)
    qU <- 1-pU
    k <- qU/pU 
    
    if (moe.sw == 1){
        rad <- e^2*(1 + k^2)^2 + k^2 * (1 - 2*e)*(1 + 2*e)
    }
    if (moe.sw == 2){
        e <- e*pU
        rad <- e^2*(1 + k^2)^2 + k^2 * (1 - 2*e)*(1 + 2*e)
    }
    x <- e*(1 + k^2) + sqrt(rad)
    x <- x / (k * (1 - 2*e))
    
    n.sam <- a * (sqrt(pU*qU)/z * log(x))^2 + 1/N
    n.sam <- 1/n.sam
    n.sam
}
