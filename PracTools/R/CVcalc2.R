CVcalc2 <- function(V=NULL, m=NULL , nbar=NULL, k=1, delta=NULL, Bsq=NULL, Wsq=NULL){
    if (any(is.null(V),is.null(m),is.null(nbar))){
         stop("V, m, and nbar must be specified.\n")
    }
    if (sum(sapply(list(Bsq, Wsq, delta), is.null)) == 3){
        stop("Either (Bsq,Wsq) or delta must be specified.\n")
    }
    if (is.null(delta)){
        if (sum(sapply(list(Bsq, Wsq), is.null)) ){
            stop("Bsq and Wsq must be specified together.\n")
        }
    }
    if (!is.null(Bsq) & is.null(Wsq) & !is.null(delta)){
        warning("Bsq specified without Wsq. delta used, Bsq ignored.\n")
    }
    if (is.null(Bsq) & !is.null(Wsq) & !is.null(delta)){
        warning("Wsq specified without Bsq. delta used, Wsq ignored.\n")
    }
    if (any(V < 0, m < 0, nbar < 0, k < 0, delta < 0, Bsq < 0, Wsq < 0)){
        stop("Illegal negative parameter specified.\n")
    }
    if (sum(sapply(list(Bsq, Wsq, delta), is.null)) == 0){
        if (delta != Bsq/(Bsq + Wsq)){
            stop("Bsq, Wsq, and delta are inconsistent. Specify only (Bsq,Wsq) or delta.\n")
        }
    }
    if (sum(sapply(list(k, Bsq, Wsq, V), is.null)) == 0){
        if (k != (Bsq + Wsq)/V){
            stop("V, k, Bsq, and Wsq are inconsistent. Specify only (V,k) or (V,Bsq,Wsq).\n")
        }
    }
    if (sum(sapply(list(Bsq, Wsq), is.null)) == 0){
        delta <- Bsq / (Bsq + Wsq)
        cv <- V/(m*nbar) * k * (1 + delta * (nbar-1))
        cv <- sqrt(cv)
    }
    if (!is.null(k)){
        cv <- V/(m*nbar) * k * (1 + delta * (nbar-1))
        cv <- sqrt(cv)
    }
    if (is.null(k)){
        k <- (Bsq + Wsq)/V
        cv <- V/(m*nbar) * k * (1 + delta * (nbar-1))
        cv <- sqrt(cv)
    }
    if (is.null(cv)) stop("Parameter combination is wrong. Check inputs.\n")
    else cv
}
