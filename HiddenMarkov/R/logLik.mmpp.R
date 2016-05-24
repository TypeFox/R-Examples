"logLik.mmpp" <-
function(object, fortran=TRUE, ...){
    tau <- object$tau[-1] - object$tau[-length(object$tau)]
    return(forwardback.mmpp(tau, object$Q, object$delta, object$lambda,
                            fortran=fortran, fwd.only=TRUE)$LL)
}


