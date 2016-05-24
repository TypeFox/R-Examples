
print.cox.ridge <-
  function (x, ...) 
  {
    cat("call:")
    cat("\n")
    print(x$call)
    cat("\n")
    nm <- unlist(attr(x$X,"dimnames")[2])
    coef <- as.vector(x$coef)
    tmp <- cbind( round(coef,4), round(exp(coef),4))
    dimnames(tmp) <- list(nm, c("coef", "exp(coef)"))
    prmatrix(tmp)
    cat("\n")
    cat("penalized log-likelihood= ", x$loglik)
    cat("\n")
    cat("Optimal penalty weight= ", x$lambda, "(converged in",x$inter.it,"internal iterations)")
    cat("\n")
    cat("Algorithm converged in", x$it, "iterations")}
