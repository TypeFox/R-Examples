nContMoe <- function(moe.sw, e, alpha=0.05, CVpop=NULL, S2=NULL, ybarU=NULL, N=Inf){
    n.sam <- NULL

    if (!(moe.sw==1 || moe.sw==2))
        stop("moe.sw must equal 1 or 2.\n")
    if (alpha <= 0 || alpha >= 1)
        stop("alpha must be in (0,1).\n")
    if (sum(sapply(list(e, N), is.null) != 0))
        stop("e and N cannot be NULL.\n")
    if (any(e <= 0) || any(e >= 1)) stop("e must be in (0,1).\n")
    if (length(alpha) > 1) stop("alpha must be scalar.\n")

    if (any(N <= 0, S2 <= 0, CVpop <= 0))
        stop("None of N, S2, or CVpop can be <= 0\n")

    if (moe.sw==1){
        if (is.null(S2))
            stop("If moe.sw=1, S2 cannot be NULL.\n")
        if (!is.null(S2) & !is.null(ybarU)){
            warning("If moe.sw=1, ybarU is ignored.\n")
        }
    }
    if (moe.sw==2){
        if (is.null(CVpop)){
            if ((is.null(S2) || is.null(ybarU))){
                stop("If moe.sw=2 and CVpop=NULL, then S2 and ybarU must be non-NULL.\n")
            }
        }
        if (!is.null(CVpop) & (!is.null(S2) || !is.null(ybarU))){
            warning("If moe.sw=2 and CVpop is non-NULL, S2 and ybarU are ignored.\n")
        }
    }
    if (sum(sapply(list(S2, ybarU, CVpop), is.null)) == 0){
        cat("S2, ybarU, and CVpop all specified. CVpop ignored.\n")
    }

    z <- qnorm(1 - alpha/2)

    if (!is.null(S2) & !is.null(ybarU)){
        CVpop <- sqrt(S2)/ybarU
    }

    if (moe.sw==1){
        n.sam <- z^2 * S2 / (e^2 + z^2*S2/N)
    }

    if (moe.sw==2){
        n.sam <- z^2 * CVpop^2 / (e^2 + z^2*CVpop^2/N)
    }

    if (is.null(n.sam)) stop("Parameter combination is wrong. Check inputs.\n")
    else n.sam
}
