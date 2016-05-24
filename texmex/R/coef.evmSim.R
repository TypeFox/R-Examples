coef.evmSim <- function(object, ...){
    res <- apply(object$param, 2, mean)
#    names(res) <- c(paste("phi:", colnames(object$X.phi)),
#        paste("xi:", colnames(object$X.xi)))
    names(res) <- names(object$map$coefficients)
    res
}
coefficients.evmSim <- coef.evmSim

