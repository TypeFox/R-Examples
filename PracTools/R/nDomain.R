nDomain <- function(CV0d=NULL, V0d=NULL, S2d=NULL, ybarUd=NULL, N=Inf, CVpopd=NULL, Pd, est.type){
    n.sam <- NULL
    if (sum(sapply(list(CV0d, V0d, S2d, ybarUd, CVpopd, Pd), is.null)) == 6){
        stop("No parameters specified\n")
    }
    if (!(est.type %in% c("total","mean"))){
        stop("est.type must be \"total\" or \"mean\".\n")
    }
    if (sum(sapply(list(CV0d, V0d), is.null)) != 1)
        stop("Only one of CV0d and V0d must be specified.\n")

    if (any(N <= 0, S2d <= 0, CV0d <= 0, V0d <=0, CVpopd <= 0))
        stop("N, S2d, CV0d, V0d, CVpopd cannot be <= 0.\n")
    if ((Pd <= 0) || (Pd > 1))
        stop("Pd must be >= 0 and < 1.\n")

    if(!is.null(V0d) & !is.null(ybarUd)){
        CV0d <- sqrt(V0d/ybarUd^2)
    }
    if (sum(sapply(list(S2d, ybarUd, CVpopd), is.null)) == 0){
        if (S2d/ybarUd^2 != CVpopd)
            stop("S2d, ybarUd, and CVpopd inconsistent. Specify only CVpopd or S2d and ybarUd.\n")
    }

    if(est.type == "total"){
        if (!is.null(CVpopd)){
           if (sum(sapply(list(CVpopd, N, CV0d, Pd), is.null)) == 0){
               Qd <- 1 - Pd
               n.sam <- (CVpopd^2 + Qd) / (Pd*CV0d^2 + (CVpopd^2 + Qd)/N)
            }
        }
        else {
            CVpopd <- sqrt(S2d / ybarUd^2)
            Qd <- 1 - Pd
            n.sam <- (CVpopd^2 + Qd) / (Pd*CV0d^2 + (CVpopd^2 + Qd)/N)
        }

    }
    if(est.type == "mean"){
        if (!is.null(CVpopd)){
           if (sum(sapply(list(CVpopd, N, CV0d, Pd), is.null)) == 0){
               n.sam <- CVpopd^2 / (Pd*CV0d^2 + CVpopd^2/N)
            }
        }
        else {
            CVpopd <- sqrt(S2d / ybarUd^2)
            n.sam <- CVpopd^2 / (Pd*CV0d^2 + CVpopd^2/N)
        }
    }
    if (is.null(n.sam)) stop("Parameter combination is wrong. Check inputs.\n")
    else n.sam
}
