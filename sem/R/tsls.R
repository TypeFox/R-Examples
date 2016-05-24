# Two-Stage Least Squares
#   John Fox

# last modified 2013-09-10 by J. Fox

tsls <- function (y, ...) {
	UseMethod("tsls")
}

tsls.default <- function (y, X, Z, w, names = NULL, ...) {
	if (is.null(w)) w <- 1
	if (any(w < 0 | is.na(w))) 
		stop("missing or negative weights not allowed")
	n <- length(y)
	p <- ncol(X)
	sqrt.w <- sqrt(w)
	invZtZ <- solve(crossprod(Z*sqrt.w))
	XtZ <- crossprod(X*w, Z)
	V <- chol2inv(chol(XtZ %*% invZtZ %*% t(XtZ)))
	b <- V %*% XtZ %*% invZtZ %*% crossprod(Z*w, y)
	residuals <- (y - X %*% b)*sqrt.w
	s2 <- sum(residuals^2)/(n - p)
	V <- s2 * V
	result <- list()
	result$n <- n
	result$p <- p
	b <- as.vector(b)
	names(b) <- names
	result$coefficients <- b
	rownames(V) <- colnames(V) <- names
	result$V <- V
	result$s <- sqrt(s2)
	result$residuals <- as.vector(residuals)
	result$response <- y
	result$model.matrix <- X
	result$instruments <- Z
	result$weights <- w
	result
}

tsls.formula <- function (formula, instruments, data, subset, weights, na.action, contrasts = NULL, ...) {
	if (missing(na.action)) 
		na.action <- options()$na.action
	m <- match.call(expand.dots = FALSE)
	if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
		m$data <- as.data.frame(data)
	response.name <- deparse(formula[[2]])
	form <- as.formula(paste(paste(response.name, collapse = ""), 
					"~", paste(deparse(formula[[3]]), collapse = ""), "+", 
					paste(deparse(instruments[[2]]), collapse = "")))
	m$formula <- form
	m$instruments <- m$contrasts <- NULL
	m[[1]] <- as.name("model.frame")
	mf <- eval(m, sys.frame(sys.parent()))
	na.act <- attr(mf, "na.action")
	w <- as.vector(model.weights(mf))
	wt.var <- if(!is.null(w)) deparse(substitute(weights)) else NULL
	Z <- model.matrix(instruments, data = mf, contrasts)
	y <- mf[, response.name]
	X <- model.matrix(formula, data = mf, contrasts)
	result <- tsls(y, X, Z, w, colnames(X))
	result$response.name <- response.name
	result$formula <- formula
	result$instruments <- instruments
	result$wt.var <- wt.var
	if (!is.null(na.act)) 
		result$na.action <- na.act
	class(result) <- "tsls"
	result
}

print.tsls <- function (x, ...) {
	cat("\nModel Formula: ")
	print(x$formula)
	cat("\nInstruments: ")
	print(x$instruments)
	if (!is.null(x$wt.var)){
		cat("\nWeights: ", x$wt.var, "\n")
	}
	cat("\nCoefficients:\n")
	print(x$coefficients)
	cat("\n")
	invisible(x)
}

summary.tsls <- function (object, digits=getOption("digits"), ...) {
    save.digits <- options(digits = digits)
    on.exit(options(save.digits))
    df <- object$n - object$p
    std.errors <- sqrt(diag(object$V))
    b <- object$coefficients
    t <- b/std.errors
    p <- 2 * (1 - pt(abs(t), df))
    table <- cbind(b, std.errors, t, p)
    rownames(table) <- names(b)
    colnames(table) <- c("Estimate", "Std. Error", "t value", 
        "Pr(>|t|)")
    
    result <- list(formula=object$formula, instruments=object$instruments, wt.var=object$wt.var,
        residuals=summary(residuals(object)), coefficients=table, digits=digits, s=object$s, df=df)
    class(result) <- "summary.tsls"
    result
}

print.summary.tsls <- function(x, ...){
    cat("\n 2SLS Estimates\n")
    cat("\nModel Formula: ")
    print(x$formula)
    cat("\nInstruments: ")
    print(x$instruments)
    if (!is.null(x$wt.var)){
        cat("\nWeights: ", x$wt.var, "\n")
    }
    cat("\nResiduals:\n")
    print(x$residuals)
    cat("\n")
    printCoefmat(x$coefficients, digits=x$digits)
    cat(paste("\nResidual standard error:", round(x$s, x$digits), 
        "on", x$df, "degrees of freedom\n\n"))
    invisible(x)
}

residuals.tsls <- function(object, ...){
    res <- object$residuals
    if (is.null(object$na.action)) 
        res
    else naresid(object$na.action, res)
    }

coef.tsls <- function(object, ...){
    object$coefficients
    }
    
fitted.tsls <- function(object, ...){
    yhat <- as.vector(object$model.matrix %*% object$coefficients)
    if (is.null(object$na.action)) 
        yhat
    else napredict(object$na.action, yhat)
    }
    
vcov.tsls <- function(object, ...) object$V
    
anova.tsls <- function(object, model.2, s2, dfe, ...){
    if(class(model.2) != "tsls") stop('requires two models of class tsls')
    s2.1 <- object$s^2
    n.1 <- object$n 
    p.1 <- object$p
    dfe.1 <- n.1 - p.1
    s2.2 <- model.2$s^2
    n.2 <- model.2$n
    p.2 <- model.2$p
    dfe.2 <- n.2 - p.2
    SS.1 <- s2.1 * dfe.1
    SS.2 <- s2.2 * dfe.2
    SS <- abs(SS.1 - SS.2)
    Df <- abs(dfe.2 - dfe.1)
    if (missing(s2)){
        s2 <- if (dfe.1 > dfe.2) s2.2 else s2.1
        f <- (SS/Df) / s2
        RSS <- c(SS.1, SS.2)
        Res.Df <- c(dfe.1, dfe.2)
        SS <- c(NA, SS)
        P <- c(NA, 1 - pf(f, Df, min(dfe.1, dfe.2)))
        Df <- c(NA, Df)
        f <- c(NA, f)
        rows <- c("Model 1", "Model 2")
        }
    else{
        f <- (SS/Df) / s2
        RSS <- c(SS.1, SS.2, s2*dfe)
        Res.Df <- c(dfe.1, dfe.2, dfe)
        SS <- c(NA, SS, NA)
        P <- c(NA, 1 - pf(f, Df, dfe), NA)
        Df <- c(NA, Df, NA)
        f <- c(NA, f, NA)
        rows <- c("Model 1", "Model 2", "Error")
        }
    table <- data.frame(Res.Df, RSS, Df, SS, f, P)
    head.1 <- paste("Model 1: ",format(object$formula), "  Instruments:", 
        format(object$instruments))
    head.2 <- paste("Model 2: ",format(model.2$formula), "  Instruments:", 
        format(model.2$instruments))
    names(table) <- c("Res.Df", "RSS", "Df", "Sum of Sq", "F", "Pr(>F)")
    row.names(table) <- rows
    structure(table, heading = c("Analysis of Variance", "", head.1, head.2, ""), 
        class = c("anova", "data.frame"))
    }
