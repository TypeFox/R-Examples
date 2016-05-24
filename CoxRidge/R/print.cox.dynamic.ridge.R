print.cox.dynamic.ridge <-
function (x, ...) 
{
    cat("call:")
    cat("\n")
    print(x$call)
    cat("\n")
    nm <- unlist(attr(x$X,"dimnames")[2])
    if(is.null(nm)){nm <- paste("Covar",1:ncol(x$X),sep="")}
    cfnames <- c( nm, paste(rep(nm, 2 - 1), rep(paste("f", 1:(ncol(x$Ft) - 1), "(t)", sep = ""), each = ncol(x$X)),  sep = ":"))
    coef <- as.vector(x$theta)
    if(!is.null(x$R)){
        coef <- c(coef,x$fixed.coef)
        fnam <- unlist(attr(x$R,"dimnames")[2])
        if(is.null(fnam)){fnam <- paste("fixed",1:ncol(x$R),sep="")}
        cfnames <- c(cfnames,fnam)}
    tmp <- cbind( round(coef,4), round(exp(coef),4))
    dimnames(tmp) <- list(cfnames, c("coef", "exp(coef)"))
    prmatrix(tmp)
    cat("\n")
    cat("penalized log-likelihood= ", x$loglik)
    cat("\n")
    cat("Optimal penalty weight= ", x$lambda, "(converged in",x$inter.it,"internal iterations)")
    cat("\n")
    cat("Algorithm converged in", x$it, "iterations")
}
