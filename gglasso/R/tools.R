coef.gglasso <- function(object, s = NULL, ...) {
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- rbind2(b0, object$beta)
    if (!is.null(s)) {
        vnames <- dimnames(nbeta)[[1]]
        dimnames(nbeta) <- list(NULL, NULL)
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)
        if(length(s) == 1)
		{
			nbeta = nbeta[, lamlist$left, drop=FALSE] * lamlist$frac +
			nbeta[, lamlist$right, drop=FALSE] * (1 - lamlist$frac)
		} else
		{
			nbeta = nbeta[, lamlist$left, drop=FALSE] %*% diag(lamlist$frac) +
			nbeta[, lamlist$right, drop=FALSE] %*% diag(1 - lamlist$frac)
		}
        dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    return(nbeta)
}


predict.gglasso <- function(object, newx, s = NULL, type = c("class", 
    "link"), ...) {
    type <- match.arg(type)
    loss <- class(object)[[2]]
    b0 <- t(as.matrix(object$b0))
    rownames(b0) <- "(Intercept)"
    nbeta <- rbind2(b0, object$beta)
    if (!is.null(s)) {
        vnames <- dimnames(nbeta)[[1]]
        dimnames(nbeta) <- list(NULL, NULL)
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)
        if(length(s) == 1)
		{
			nbeta = nbeta[, lamlist$left, drop=FALSE] * lamlist$frac +
			nbeta[, lamlist$right, drop=FALSE] * (1 - lamlist$frac)
		} else
		{
			nbeta = nbeta[, lamlist$left, drop=FALSE] %*% diag(lamlist$frac) +
			nbeta[, lamlist$right, drop=FALSE] %*% diag(1 - lamlist$frac)
		}
        dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
    }
    if (is.null(dim(newx))) newx = matrix(newx, nrow = 1)
    nfit <- as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
    if (loss %in% c("logit", "sqsvm", "hsvm")) {
        switch(type, link = nfit, class = ifelse(nfit > 0, 1, -1))
    } else {
        nfit
    }
}


print.gglasso <- function(x, digits = max(3, getOption("digits") - 
    3), ...) {
    cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Df = x$df, Lambda = signif(x$lambda, digits)))
}


coef.cv.gglasso <- function(object, s = c("lambda.1se", "lambda.min"), 
    ...) {
    if (is.numeric(s)) 
        lambda <- s else if (is.character(s)) {
        s <- match.arg(s)
        lambda <- object[[s]]
    } else stop("Invalid form for s")
    coef(object$gglasso.fit, s = lambda, ...)
}

predict.cv.gglasso <- function(object, newx, s = c("lambda.1se", 
    "lambda.min"), ...) {
    if (is.numeric(s)) 
        lambda <- s else if (is.character(s)) {
        s <- match.arg(s)
        lambda <- object[[s]]
    } else stop("Invalid form for s")
    predict(object$gglasso.fit, newx, s = lambda, ...)
} 
