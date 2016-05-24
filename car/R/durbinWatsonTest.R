#-------------------------------------------------------------------------------
# Revision history:
# 2009-09-28 by J. Fox (renamed)
#-------------------------------------------------------------------------------

# generalized Durbin-Watson statistic (J. Fox)

durbinWatsonTest <- function(model, ...){
	UseMethod("durbinWatsonTest")
}

durbinWatsonTest.lm <- function(model, max.lag=1, simulate=TRUE, reps=1000, 
	method=c("resample","normal"), 
	alternative=c("two.sided", "positive", "negative"), ...){
	method <- match.arg(method)
	alternative <- if (max.lag == 1) match.arg(alternative)
		else "two.sided"
	residuals <- residuals(model)
	if (any(is.na(residuals))) stop ('residuals include missing values')
	n <- length(residuals)
	r <- dw <-rep(0, max.lag)
	den <- sum(residuals^2)
	for (lag in 1:max.lag){
		dw[lag] <- (sum((residuals[(lag+1):n] - residuals[1:(n-lag)])^2))/den
		r[lag] <- (sum(residuals[(lag+1):n]*residuals[1:(n-lag)]))/den
	}
	if (!simulate){
		result <- list(r=r, dw=dw)
		class(result) <- "durbinWatsonTest"
		result
	}
	else {
		S <- summary(model)$sigma
		X <- model.matrix(model)
		mu <- fitted.values(model)
		Y <- if (method == "resample") 
				matrix(sample(residuals, n*reps, replace=TRUE), n, reps) + matrix(mu, n, reps)
			else matrix(rnorm(n*reps, 0, S), n, reps) + matrix(mu, n, reps)
		E <- residuals(lm(Y ~ X - 1))
		DW <- apply(E, 2, durbinWatsonTest, max.lag=max.lag)
		if (max.lag == 1) DW <- rbind(DW)
		p <- rep(0, max.lag)
		if (alternative == 'two.sided'){
			for (lag in 1:max.lag) {
				p[lag] <- (sum(dw[lag] < DW[lag,]))/reps
				p[lag] <- 2*(min(p[lag], 1 - p[lag]))
			}
		}
		else if (alternative == 'positive'){
			for (lag in 1:max.lag) {
				p[lag] <- (sum(dw[lag] > DW[lag,]))/reps
			}
		}
		else {
			for (lag in 1:max.lag) {
				p[lag] <- (sum(dw[lag] < DW[lag,]))/reps
			}
		}
		result <- list(r=r, dw=dw, p=p, alternative=alternative)
		class(result)<-"durbinWatsonTest"
		result
	}
}

durbinWatsonTest.default <- function(model, max.lag=1, ...){
	# in this case, "model" is the residual vectors
	if ((!is.vector(model)) || (!is.numeric(model)) ) stop("requires vector of residuals")
	if (any(is.na(model))) stop ('residuals include missing values')
	n <-  length(model)
	dw <- rep(0, max.lag)
	den <- sum(model^2)
	for (lag in 1:max.lag){
		dw[lag] <- (sum((model[(lag+1):n] - model[1:(n-lag)])^2))/den
	}
	dw
}

print.durbinWatsonTest <- function(x, ...){
	max.lag <- length(x$dw)
	result <- if (is.null(x$p)) cbind(lag=1:max.lag,Autocorrelation=x$r, "D-W Statistic"=x$dw)
		else cbind(lag=1:max.lag,Autocorrelation = x$r, "D-W Statistic" = x$dw, 
				"p-value"= x$p)
	rownames(result) <- rep("", max.lag)
	print(result)
	cat(paste(" Alternative hypothesis: rho", if(max.lag > 1) "[lag]" else "",
			c(" != ", " > ", " < ")[which(x$alternative == c("two.sided", "positive", "negative"))],
			"0\n", sep=""))
	invisible(x)
}

dwt <- function(...) durbinWatsonTest(...)
