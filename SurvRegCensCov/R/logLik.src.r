logLik.src <- function(object, ...){
     p1 <- paste("'log Lik.' ", round(object$loglik, 5), " (df=", object$d, ")", sep = "")
     cat(p1)
     invisible(object)
}