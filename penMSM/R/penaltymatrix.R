penaltymatrix <- function(lambda, PSM, beta, w, constant){
    if(length(lambda) != nrow(PSM)){
        stop("check dimensions of lambda and penalty-structure-matrix.")
    }
    p <- ncol(PSM)
    L <- nrow(PSM)
    plambda <- matrix(nrow=p, ncol=p, 0)
    for(l in 1:L){
        helpobject <- plmatrix(psv = PSM[l, ], beta = beta, constant = constant)
        plambda <- plambda + lambda[l] * w[l] * helpobject
    }
    return(plambda)
}
