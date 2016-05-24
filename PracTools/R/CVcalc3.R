CVcalc3 <- function(V=NULL, m=NULL , nbar=NULL, qbar=NULL, k1=1, k2=1, 
                delta1=NULL, delta2=NULL, Bsq=NULL, Wsq=NULL, W2sq=NULL, W3sq=NULL){
    if (any(is.null(V),is.null(m),is.null(nbar),is.null(qbar))){
         stop("V, m, nbar, qbar must be specified.\n")
    }
    if (sum(sapply(list(Bsq, Wsq, delta1), is.null)) == 3){
        stop("Either (Bsq,Wsq) or delta1 must be specified.\n")
    }
    if (sum(sapply(list(W2sq, W3sq, delta2), is.null)) == 3){
        stop("Either (W2sq,W3sq) or delta2 must be specified.\n")
    }
    if (is.null(delta1)){
        if (sum(sapply(list(Bsq, Wsq), is.null)) ){
            stop("Bsq and Wsq must be specified together.\n")
        }
    }
    if (is.null(delta2)){
        if (sum(sapply(list(W2sq, W3sq), is.null)) ){
            stop("W2sq and W3sq must be specified together.\n")
        }
    }

    if (!is.null(Bsq) & is.null(Wsq) & !is.null(delta1)){
        warning("Bsq specified without Wsq. delta1 used, Bsq ignored.\n")
    }
    if (is.null(Bsq) & !is.null(Wsq) & !is.null(delta1)){
        warning("Wsq specified without Bsq. delta1 used, Wsq ignored.\n")
    }
    if (!is.null(W2sq) & is.null(W3sq) & !is.null(delta2)){
        warning("W2sq specified without W3sq. delta2 used, W2sq ignored.\n")
    }
    if (is.null(W2sq) & !is.null(W3sq) & !is.null(delta2)){
        warning("W3sq specified without W2sq. delta2 used, W3sq ignored.\n")
    }

    if (any(V < 0, m < 0, nbar < 0, qbar < 0, k1 < 0, k2 < 0, 
            delta1 < 0, delta2 < 0, Bsq < 0, Wsq < 0, W2sq < 0, W3sq < 0)){
        stop("Illegal negative parameter specified.\n")
    }
    if (sum(sapply(list(Bsq, Wsq, delta1), is.null)) == 0){
        if (delta1 != Bsq/(Bsq + Wsq)){
            stop("Bsq, Wsq, and delta1 are inconsistent. Specify only (Bsq,Wsq) or delta1.\n")
        }
    }
    if (sum(sapply(list(W2sq, W3sq, delta2), is.null)) == 0){
        if (delta2 != W2sq/(W2sq + W3sq)){
            stop("W2sq, W3sq, and delta2 are inconsistent. Specify only (W2sq,W3sq) or delta2.\n")
        }
    }

    if (sum(sapply(list(k1, Bsq, Wsq, V), is.null)) == 0){
        if (k1 != (Bsq + Wsq)/V){
            stop("V, k1, Bsq, and Wsq are inconsistent. Specify only (V,k) or (V,Bsq,Wsq).\n")
        }
    }
    
    if (sum(sapply(list(k2, Bsq, Wsq, V), is.null)) == 0){
        if (k2 != (W2sq + W3sq)/V){
            stop("V, k2, W2sq, and W3sq are inconsistent. Specify only (V,k2) or (V,W2sq,W3sq).\n")
        }
    }
    if (sum(sapply(list(Bsq, Wsq), is.null)) == 0){
        delta1 <- Bsq / (Bsq + Wsq)
        if  (sum(sapply(list(W2sq, W3sq), is.null)) == 0) {delta2 <- W2sq / (W2sq + W3sq)}
        if (is.null(k1)) {k1 <- (Bsq + Wsq)/V}
        if (is.null(k2)) {k2 <- (W2sq + W3sq)/V}
        cv <- V/(m*nbar*qbar) * ( k1*delta1*nbar*qbar + k2*(1 + delta2 * (qbar-1)) )
        cv <- sqrt(cv)
    }
    if (!is.null(k1) & !is.null(k2)){
        if (is.null(delta1)) {delta1 <- Bsq / (Bsq + Wsq)}
        if (is.null(delta2)) {delta2 <- W2sq / (W2sq + W3sq)}
        cv <- V/(m*nbar*qbar) * ( k1*delta1*nbar*qbar + k2*(1 + delta2 * (qbar-1)) )
        cv <- sqrt(cv)
    }
    if (is.null(cv)) stop("Parameter combination is wrong. Check inputs.\n")
    else cv
}
