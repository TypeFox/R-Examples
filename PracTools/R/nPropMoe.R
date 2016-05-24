
nPropMoe <- function(moe.sw, e, alpha=0.05, pU, N=Inf){
    n.sam <- NULL

    if (!(moe.sw==1 || moe.sw==2))
        stop("moe.sw must equal 1 or 2.\n")
    if (alpha <= 0 || alpha >= 1)
        stop("alpha must be in (0,1).\n")
    if (sum(sapply(list(e, N, pU), is.null) != 0))
        stop("e, N, and pU cannot be NULL.\n")
    if (any(e <= 0) || any(e >= 1)) stop("e must be in (0,1).\n")
    if (any(pU <= 0) || any(pU >= 1)) stop("pU must be in (0,1).\n")
    if (N <= 0) stop("N must be positve.\n")
    if (length(alpha) > 1) stop("alpha must be scalar.\n")

    if (N == Inf) {a <- 1}
            else {a <- N/(N-1)}

    z <- qnorm(1 - alpha/2)
    qU <- 1-pU

    if (moe.sw==1){
        n.sam <- a * z^2 *pU*qU / (e^2 + z^2*pU*qU/(N-1) )
    }

    if (moe.sw==2){
        n.sam <- a * z^2 * qU/pU / (e^2 + z^2*qU/pU/(N-1) )
    }

    if (is.null(n.sam)) stop("Parameter combination is wrong. Check inputs.\n")
    else n.sam
}
