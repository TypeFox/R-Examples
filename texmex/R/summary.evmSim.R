summary.evmSim <- function(object, ...){
   co <- coef(object)
   se <- apply(object$param, 2, sd)
   res <- cbind(co, se)
   dimnames(res) <- list(names(co), c("Posterior mean", "SD"))
   res <- list(object$map$family, res)
   oldClass(res) <- "summary.evm.sim"
   res
}

print.summary.evmSim <- function(x, ...){
   print(x[[1]], verbose=FALSE, ...)
   cat("\nPosterior summary:\n")
   print(unclass(x[[2]]))
   invisible()
}