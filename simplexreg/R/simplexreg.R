simplexreg <-
function(formula, data, subset, na.action, type = c("homo", "hetero"), 
	link = c("logit", "probit", "cloglog", "neglog"), corr = "Ind", id = NULL, 
	control = simplexreg.control(...), model = TRUE, y = TRUE, x = TRUE, ...)
{   
   	call <- match.call()
   	if (missing(data))
   	   	data <- environment(formula)
   	if (missing(type))
   	   	type <- "homo"
   	if (missing(link))
   	   	link <- "logit"	
   	mf <- match.call(expand.dots = FALSE)
   	m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
   	mf <- mf[c(1L, m)]
   	mf$drop.unused.levels <- TRUE
   	oformula <- as.formula(formula)
   	options(warn = -1)
   	formula <- as.Formula(formula)
   	if (length(formula)[2L] < 2L){
   	   	formula <- as.Formula(formula(formula), ~1)	
   	}
   	else {
   	   	if (length(formula)[2L] == 2L) 
   	   	   	formula <- Formula(formula(formula, rhs = 1:2))
   	   	else{
   	   	   	formula <- Formula(formula(formula, rhs = 1:3))
		}
   	}
   	mf$formula <- formula
   	mf[[1L]] <- as.name("model.frame")
   	mf <- eval(mf, parent.frame())
   	mt <- terms(formula, data = data)
   	mtX <- terms(formula, data = data, rhs = 1L)
   	mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
   	mtT <- delete.response(terms(formula, data = data, rhs = 3L))
   	Y <- model.response(mf, "numeric")
   	X <- model.matrix(mtX, mf)
   	Z <- model.matrix(mtZ, mf)
   	T <- model.matrix(mtT, mf)
   	T <- as.vector(T[,-1])
   	options(warn = 0)
   	if (!missing(id)){
   	   	namevar <- deparse(substitute(id))
   	   	if (!missing(data))
   	   	   	id <- data[, namevar]
   	   	if (!missing(subset))	
   	   	   	id <- as.vector(id[subset])	
   	   	id <- as.factor(id)
   	   	if (length(id) != length(Y))
   	   	   	stop("id must have the same length as Y")
   	}	
   	names = list(x = dimnames(X)[[2]], z = dimnames(Z)[[2]])
   	if (length(Y) < 1)
   	   	stop("empty model")
   	if (min(Y) <= 0 || max(Y) >= 1)
   	   	stop("observations must be in (0, 1)")
   	result <- simplexreg.fit(y = Y, x = X, z = Z, t = T, type = type, link = link, corr = corr, 
   	   	id = id, control = control)
	if (is.null(result)){
		warning("Independent correlation structure is used in the marginal model")
		corr <- "Ind"
		result <- simplexreg.fit(y = Y, x = X, z = Z, t = T, type = type, link = link, corr = corr, 
			id = id, control = control)
	}
   	result$terms <- list(mean = mtX, dispersion = mtZ)
   	result$levels <- list(mean = .getXlevels(mtX, mf), dispersion = .getXlevels(mtZ, mf))
   	result$contrasts <- list(mean = attr(X, "contrasts"), dispersion = attr(Z, "contrasts"))
   	dimnames(result$fixef)[[1]] <- names$x
   	if (type == "hetero")
   	   	dimnames(result$dispar)[[1]] <- names$z
   	meanmu <- result$meanmu
   	disper <- result$Disper
   	result$call <- call
   	result$formula <- oformula
   	result$link <- link
   	result$type <- type
   	result$n <- length(Y)
   	result$corr <- corr
   	if (model)
   	   	result$model <- mf
   	if (y)
   	   	result$y <- Y
   	if (x)
		result$x <- list(mean = X, dispersion = Z, time = T, id = id)
   	class(result) <- "simplexreg"
   	return(result)
}

simplexreg.control <- function(maxit = 200, beta = NULL, gamma = NULL, alpha = NULL, 
   	tol = 1e-6, ...)
{
   	rval <- list(maxit = maxit, beta = beta, gamma = gamma, alpha = alpha, tol = tol)
   	rval
}

print.simplexreg <- 
function(x, digits = max(3, getOption("digits") - 3), ...)
{
   	cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

    if(length(x$fixef[,1])) {
   	   	cat(paste("Coefficients (mean model with ", x$link, " link):\n", sep = ""))
   	   	print.default(format(x$fixef[,1], digits = digits), print.gap = 2, quote = FALSE)
   	   	cat("\n")
   	} 
   	else cat("No coefficients (in mean model)\n\n")
   	if(x$type != "hetero") {
   	   	cat(paste("Dispersion:\n", sep = ""))
   	   	names(x$Disper) <- "(sigma^2)"
   	   	print.default(format(x$Disper, digits = digits), print.gap = 2, quote = FALSE)
   	   	cat("\n")
   	}
   	else{
   	   	cat(paste("Coefficients (dispersion model with log link):\n", sep = ""))
   	   	print.default(format(x$dispar[,1], digits = digits), print.gap = 2, quote = FALSE)
   	   	cat("\n")
   	}
	
   	if(x$corr != "Ind"){
   	   	cat(paste("Coefficients (the correlation):\n", sep = ""))
   	   	names(x$autocor) <- rep("(rho)", 4)
   	   	print.default(format(x$autocor[3], digits = digits), print.gap = 2, quote = FALSE)
   	   	cat("\n")
   	}
   	invisible(x)
}

summary.simplexreg <- 
function(object, type = "stdPerr", ...)
{
   	corr <- object$corr

   	## residuals
   	type <- match.arg(type, c("appstdPerr", "stdPerr", "stdscor"))
   	if (type == "appstdPerr")
   	   	object$residuals <- object$appstdPerr
   	else if (type == "stdPerr")
   	   	object$residuals <- object$stdPerr
   	else 
   	   	object$residuals <- object$stdscor
   	object$residuals.type <- type

   	## extend coefficient table
   	type <- object$type
   	if (type == "homo"){
   	   	cf <- as.vector(object$fixef[,1])
   	   	se <- as.vector(object$fixef[,2])
   	   	k <- length(object$fixef[,1])
   	}
   	else{
   	   	cf <- as.vector(c(object$fixef[,1], object$dispar[,1]))
   	   	se <- as.vector(c(object$fixef[,2], object$dispar[,2]))
   	   	k <- length(object$fixef[,1])
   	   	m <- length(object$dispar[,1])
   	}
   	if (corr != "Ind"){
   	   	cf <- c(cf, object$autocor[c(1, 3)])
   	   	se <- c(se, object$autocor[c(2, 4)])
   	}
   	cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
   	colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
   	if (type == "homo"){
   	   	if (corr != "Ind"){
   	   	   	cf <- list(mean = cf[seq.int(length.out = k), , drop = FALSE], correlation = 
   	   	   	   	cf[seq.int(length.out = 2) + k, , drop = FALSE])
   	   	   	rownames(cf$mean) <- names(object$fixef[,1])
   	   	   	rownames(cf$correlation) <- c("alpha", "rho")
   	   	}
   	   	else{
   	   	   	cf <- list(mean = cf)
   	   	   	rownames(cf$mean) <- names(object$fixef[,1])
   	   	}
   	}
   	else{
   	   	if (corr != "Ind"){	
   	   	   	cf <- list(mean = cf[seq.int(length.out = k), , drop = FALSE], dispersion = 
   	   	   	   	cf[seq.int(length.out = m) + k, , drop = FALSE], correlation = 
   	   	   	   	cf[seq.int(length.out = 2) + m + k, , drop = FALSE])
   	   	   	rownames(cf$mean) <- names(object$fixef[,1])			
   	   	   	if (type == "hetero")
   	   	   	   	rownames(cf$dispersion) <- names(object$dispar[,1])
   	   	   	else
   	   	   	   	rownames(cf$dispersion) <- "(Intercept)"
   	   	   	rownames(cf$correlation) <- c("alpha", "rho")
   	   	}
   	   	else{
   	   	   	cf <- list(mean = cf[seq.int(length.out = k), , drop = FALSE], dispersion = 
   	   	   	   	cf[seq.int(length.out = m) + k, , drop = FALSE])
   	   	   	rownames(cf$mean) <- names(object$fixef[,1])
   	   	   	if (type == "hetero")
   	   	   	   	rownames(cf$dispersion) <- names(object$dispar[,1])
   	   	   	else
   	   	   	   	rownames(cf$dispersion) <- "(Intercept)"	
   	   	}
   	}
   	object$coefficients <- cf

   	## delete some slots
   	object$predict <- object$meanmu <- object$model <- object$y <- object$terms <- object$contrasts <-
   	object$levels <- object$x <- object$stdPerr <- object$appstdPerr <- object$adjvar <- 
   	object$fixef <- object$dispar <- object$Dispersion <- object$adjvar <- NULL

   	## return
   	class(object) <- "summary.simplexreg"
   	object
}

print.summary.simplexreg <- 
function(x, digits = max(3, getOption("digits") - 3), ...)
{
   	cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

   	types <- c("appstdPerr", "stdPerr", "stdscor")
   	Types <- c("approximate Pearson residuals", "standard Pearson residuals", "standardised score residuals")
   	cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
   	print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
   	   	.Names = c("Min", "1Q", "Median", "3Q", "Max")))

   	if(NROW(x$coefficients$mean)) {
   	   	cat(paste("\nCoefficients (mean model with ", x$link, " link):\n", sep = ""))
   	   	printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
   	} 

   	if(NROW(x$coefficients$dispersion)) {
   	   	cat(paste("\nCoefficients (dispersion model with log link):\n", sep = ""))
   	   	printCoefmat(x$coefficients$dispersion, digits = digits, signif.legend = FALSE)
   	} 
	
   	if(NROW(x$coefficients$correlation)) {
   	   	cat(paste("\nCoefficients (correlation):\n", sep = ""))
   	   	printCoefmat(x$coefficients$correlation, digits = digits, signif.legend = FALSE)
   	} 

   	if (getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
   	   	cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

   	corr <- x$corr
   	type <- x$type
	
   	if (corr == "Ind"){
   	   	df <- x$n - NROW(x$coefficients$mean)
   	   	cat("\nLog-likelihood:", formatC(x$loglik, digits = digits))
   	   	cat(",  p-value:", 1 - pchisq(x$deviance, df), "\n")
   	   	cat("Deviance:", x$deviance, "\n")
   	}
   	else{
   	   	cat("\nOverall Deviance:", sum(x$deviance), "\n")
   	}
	
   	iter <- x$iter
   	cat("Number of Fisher Scoring iterations: ", iter, "\n")
   	invisible(x)
}

AIC.simplexreg <- 
function(object, ..., k = 2){
	if (object$corr != "Ind")
		stop("Calculating AIC from GEE models is not supported")
	if (object$type != "hetero")
		df <- nrow(object$fixef) + 1
	else 
		df <- nrow(object$fixef) + nrow(object$dispar)
	return(2 * object$loglike - 2 * df)
}

coef.simplexreg <-
function(object, ...){
   	type <- object$type
   	if (type == "homo"){
   	   	cf <- as.vector(object$fixef[,1])
   	   	se <- as.vector(object$fixef[,2])
   	}
   	else{
   	   	cf <- as.vector(c(object$fixef[,1], object$dispar[,1]))
   	   	se <- as.vector(c(object$fixef[,2], object$dispar[,2]))
   	}
   	corr <- object$corr
   	if (corr != "Ind"){
   	   	cf <- as.vector(c(cf, object$auto[c(1,3)]))
   	   	se <- as.vector(c(se, object$auto[c(2,4)]))
   	}
   	cf <- cbind(cf, se)
   	colnames(cf) <- c("Estimate", "Std. Error")
   	if (type == "homo")
   	   	Rname <- names(object$fixef[,1])
   	else
   	   	Rname <- c(names(object$fixef[,1]), names(object$dispar[,1]))
   	if (corr != "Ind")
   	   	Rname <- c(Rname, "alpha", "rho")
   	rownames(cf) <- Rname	
	return(cf)
}

residuals.simplexreg <- function(object, type = c("appstdPerr", "stdPerr", "stdscor", 
   	"adjvar"), ...)
{
   	## residuals
   	type <- match.arg(type, c("appstdPerr", "stdPerr", "stdscor", "adjvar"))

   	if (type == "appstdPerr")
   	   	res <- object$appstdPerr
   	else if (type == "stdPerr")
   	   	res <- object$stdPerr
   	else if (type == "stdscor")
   	   	res <- object$stdscor
   	else
   	   	res <- object$adjvar

   	return(res)
}

predict.simplexreg <- 
function(object, newdata = NULL, type = c("response", "dispersion"), na.action, ...)
{
   	type <- match.arg(type)
   	if(missing(newdata)) {
   	   	if (type == "response")
   	   	   	result <- object$meanmu
   	   	else
   	   	   	result <- object$dispersion
   	} 
   	else {
   	   	mf <- model.frame(newdata)
   	   	if (type == "response"){
   	   	   	mf <- model.frame(delete.response(object$terms[[1]]), newdata, na.action = na.action, 
   	   	   	   	xlev = object$levels[[1]])
   	   	   	X <- model.matrix(delete.response(object$terms$mean), mf, contrasts = object$contrasts$mean)
			result <- X %*% object$fixef[,1]
			if (object$link == "logit")
				result <- hf1(result)
			else if (object$link == "probit")
				result <- hf2(result)
			else if (object$link == "cloglog")
				result <- hf3(result)
			else
				result <- hf4(result)
		}
   	   	else{
   	   	   	if (object$type == "homo")
   	   	   	   	stop("Prediction for dispersion is not supported for homogeneous model")
			mf <- model.frame(delete.response(object$terms[[2]]), newdata, na.action = na.action, 
   	   	   	   	xlev = object$levels[[2]])
			Z <- model.matrix(object$terms$dispersion, mf, contrasts = object$contrasts$dispersion)
			result <- exp(Z %*% object$dispar[,1])
		}
    
    return(result)
	}
}

plot.simplexreg <- 
function(x, type = c("residuals", "corr", "GOF"), res = "adjvar", lag = 1, ...)
{
   	type = match.arg(type)
   	dots <- list(...)
   	if (type == "residuals"){
   	   	res <- match.arg(res, c("appstdPerr", "stdPerr", "stdscor", "adjvar"))
   	   	if (res == "appstdPerr"){
   	   	   	res <- x$appstdPerr
   	   	   	if (is.null(dots$xlab)){
   	   	   	   	if (is.null(dots$ylab))
   	   	   	   	   	plot(res~x$meanmu, xlab = "Estimated Mean", 
   	   	   	   	   	    ylab = "Approximate Pearson Residuals", ...)
   	   	   	   	else
   	   	   	   	   	plot(res~x$meanmu, xlab = "Estimated Mean", ...)
   	   	   	}
   	   	   	else{
   	   	   	   	if (is.null(dots$ylab))
   	   	   	   	   	plot(res~x$meanmu, ylab = "Approximate Pearson Residuals", ...)
   	   	   	   	else
   	   	   	   	   	plot(res~x$meanmu, ...)
   	   	   	}
   	   	   	abline(h = c(-1.96, 1.96))
   	   	}
   	   	else if (res == "stdPerr"){
   	   	   	res <- x$stdPerr
   	   	   	if (is.null(dots$xlab)){
   	   	   	   	if (is.null(dots$ylab))
   	   	   	   	   	plot(res~x$meanmu, xlab = "Estimated Mean", 
   	   	   	   	   	   	ylab = "Pearson Residuals", ...)
   	   	   	   	else
   	   	   	   	   	plot(res~x$meanmu, xlab = "Estimated Mean", ...)
   	   	   	}
   	   	   	else{
   	   	   	   	if (is.null(dots$ylab))
   	   	   	   	   	plot(res~x$meanmu, ylab = "Pearson Residuals", ...)
   	   	   	   	else
   	   	   	   	   	plot(res~x$meanmu, ...)
   	   	   	}
   	   	   	abline(h = c(-1.96, 1.96))
   	   	}
   	   	else if (res == "stdscor"){
   	   	   	res <- x$stdscor
   	   	   	if (is.null(dots$xlab)){
   	   	   	   	if (is.null(dots$ylab))
   	   	   	   	   	plot(res~x$meanmu, xlab = "Estimated Mean", 
   	   	   	   	   	   	ylab = "Standardised Score Residuals", ...)
   	   	   	   	else
   	   	   	   	   	plot(res~x$meanmu, xlab = "Estimated Mean", ...)
   	   	   	}
   	   	   	else{
   	   	   	   	if (is.null(dots$ylab))
   	   	   	   	   	plot(res~x$meanmu, ylab = "Approximate Pearson Residuals", ...)
   	   	   	   	else
   	   	   	   	   	plot(res~x$meanmu, ...)
   	   	   	}
   	   	   	abline(h = c(-1.96, 1.96))
   	   	}
   	   	else{
   	   	   	res <- x$adjvar 
   	   	   	if (x$link == "logit")
   	   	   	   	label <- "Logit Linear Predictor"
   	   	   	else if (x$link == "probit")
   	   	   	   	label <- "Probit Linear Predictor"
   	   	   	else if (x$link == "cloglog")
   	   	   	   	label <- "Cloglog Predictor"
   	   	   	else
   	   	   	   	label <- "Neglog Predictor"
			if (is.null(dots$xlab)){
   	   	   	   	if (is.null(dots$ylab))
   	   	   	   	   	plot(res~x$meanmu, xlab = label, ylab = "Adjusted Dependent Variables", ...)
   	   	   	   	else
   	   	   	   	   	plot(res~x$meanmu, xlab = label, ...)
   	   	   	}
   	   	   	else{
   	   	   	   	if (is.null(dots$ylab))
   	   	   	   	   	plot(res~x$meanmu, ylab = "Adjusted Dependent Variables", ...)
   	   	   	   	else
   	   	   	   	   	plot(res~x$meanmu, ...)
   	   	   	}
   	   	   	plot(res~x$predict, ...)
   	   	   	abline(a = -1.96/sqrt(2), b = 1)
   	   	   	abline(a = 1.96/sqrt(2), b = 1)
   	   	}
   	}
	
   	else if (type == "corr"){
   	   	id <- x$x$id
   	   	if (is.null(id))
   	   	   	stop("The clusters should be specified")
   	   	res <- x$stdscor
   	   	T <- x$x$time
   	   	if (is.null(T))
   	   	   	stop("The time covariate should be specified")
   	   	mm <- cormm(id, T, res)
   	   	if (is.null(dots$xlab)){
   	   	   	if (is.null(dots$ylab))
   	   	   	   	plotmm(mm, lag, xlab = "Standardised Score Residuals", 
   	   	   	   	   	ylab = "Standardised Score Residuals", ...)
   	   	   	else
   	   	   	   	plotmm(mm, lag, xlab = "Standardised Score Residuals", ...)
   	   	}
   	   	else{
   	   	   	if (is.null(dots$ylab))
   	   	   	   	plotmm(mm, lag, ylab = "Standardised Score Residuals", ...)
   	   	   	else
   	   	   	   	plotmm(mm, lag, ...)
   	   	}
   	}
	
   	else{
   	   	corr <- x$corr
   	   	if (corr == "Ind")
   	   	   	stop("Partial deviances are defined on gee models")
   	   	T <- x$x$time
   	   	if (is.null(T))
   	   	   	stop("The time covariate should be specified")
   	   	p <- dim(x$x$mean)[2]
   	   	if (is.null(p))
   	   	   	stop("The mean covariate should be specified")
   	   	time <- sort(unique(T))
   	   	tlen <- length(time)
   	   	df <- rep(0, tlen)
   	   	for (i in 1:tlen){
   	   	   	df[i] <- sum(T == time[i])
   	   	}
   	   	devi <- x$deviance
   	   	time <- time[df > p]
   	   	df <- df[df > p]
   	   	df <- df - p
   	   	tlen <- length(time)
   	   	pardevi <- rep(0, tlen)
   	   	for (i in 1:tlen){
   	   	   	pardevi[i] = sum(devi[T == time[i]])
   	   	}
   	   	if (is.null(dots$ylab))
   	   	   	plot(qchisq(0.95, df = df)~time, type = "l", lty = 3, 
   	   	   	   	ylab = "Partial Deviance", ...)
   	   	else
   	   	   	plot(qchisq(0.95, df = df)~time, type = "l", lty = 3, ...)
   	   	for (i in 1:tlen){
   	   	   	ablineclip(v = time[i], y1 = 0, y2 = pardevi[i])
   	   	}
   	}
}