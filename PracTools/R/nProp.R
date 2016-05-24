nProp <-
function(CV0=NULL, V0=NULL, pU=NULL, N=Inf){
    n.sam <- NULL
    if (sum(sapply(list(N, pU), is.null) != 0))
        stop("N and pU cannot be NULL.\n")
    if (sum(sapply(list(CV0, V0), is.null)) != 1)
        stop("Either CV0 or V0 must be specified.\n")
    if (any(pU <= 0) | any(pU >= 1)) stop("pU must be in (0,1).\n")

    if (any(N <= 0, CV0 <= 0, V0 <=0)) 
        stop("N, CV0, and V0 cannot be <= 0.\n")
    
    if (sum(sapply(list(pU, N, CV0), is.null)) == 0){
        if (N == Inf) {a <- 1}
            else {a <- N/(N-1)}
        qU <- 1-pU
        n.sam <- a * qU/pU / (CV0^2 + qU/pU/(N-1))
    }
    if (sum(sapply(list(pU, N, V0), is.null)) == 0){
        if (N == Inf) {a <- 1}
            else {a <- N/(N-1)}
        qU <- 1-pU
        n.sam <- a * pU*qU / (V0 + pU*qU/(N-1))
    }
    
    if (is.null(n.sam)) stop("Parameter combination is wrong. Check inputs.\n")
    else n.sam
}

