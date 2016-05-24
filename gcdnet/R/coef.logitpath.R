coef.logitpath <- function(object, s = NULL, type = c("coefficients", 
    "nonzero"), ...) {
    type <- match.arg(type)
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- rbind2(b0, object$beta)
    if (!is.null(s)) {
        vnames <- dimnames(nbeta)[[1]]
        dimnames(nbeta) <- list(NULL, NULL)
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)
        nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + 
            nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
        dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    if (type == "coefficients") 
        return(nbeta)
    if (type == "nonzero") 
        return(nonzero(nbeta[-1, , drop = FALSE], bystep = TRUE))
} 
