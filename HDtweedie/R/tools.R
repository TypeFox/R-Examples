coef.HDtweedie <- function(object, s = NULL, ...) {
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- rbind2(b0, object$beta)
    if (!is.null(s)) {
        vnames <- dimnames(nbeta)[[1]]
        dimnames(nbeta) <- list(NULL, NULL)
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)
        nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac + nbeta[, 
            lamlist$right, drop = FALSE] * (1 - lamlist$frac)
        dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    return(nbeta)
}


predict.HDtweedie <- function(object, newx, s = NULL, 
    type = c("response", "link"), ...) {
    if (is.vector(newx)) newx <- t(as.matrix(newx))
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
    nfit <- as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
    switch(type, link = nfit, response = exp(nfit))
} 


print.HDtweedie <- function(x, digits = max(3, getOption("digits") - 
    3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Df = x$df, Lambda = signif(x$lambda, digits)))
}


coef.cv.HDtweedie <- function(object, s = c("lambda.1se", "lambda.min"), 
    ...) {
    if (is.numeric(s)) 
        lambda <- s else if (is.character(s)) {
        s <- match.arg(s)
        lambda <- object[[s]]
    } else stop("Invalid form for s")
    coef(object$HDtweedie.fit, s = lambda, ...)
}

predict.cv.HDtweedie <- function(object, newx, s = c("lambda.1se", 
    "lambda.min"), ...) {
    if (is.numeric(s)) 
        lambda <- s else if (is.character(s)) {
        s <- match.arg(s)
        lambda <- object[[s]]
    } else stop("Invalid form for s")
    predict(object$HDtweedie.fit, newx, s = lambda, ...)
} 
