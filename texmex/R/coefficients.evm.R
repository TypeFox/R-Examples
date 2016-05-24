coefficients.evmOpt <-
function(object, ...){
    res <- object$coefficients
#    if (length(res) == 2)
#    names(res) <- c("phi", "xi")
    res
}

coef.evmOpt <- coefficients.evmOpt
#function(object, ...) {
#  coefficients.evm(object)
#}

