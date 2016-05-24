
print.summary.qme <- function(x, digits= max(3, getOption("digits") - 3), ...)
{
	cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
	if (length(x$coefficients))	 
	{
 		cat("Coefficients:\n")
		printCoefmat(x$coefficients, digits = digits, print.gap=2.5)
 	}
	else cat("No coefficients\n")
	cat("\n")
	if(length(x$counts))
	{
		cat("Number of iterations:\n")
		print(x$counts)
	}
  cat("\n")
  if(length(x$cval))
  {
    cat("Threshold information:\n")
    print.data.frame(x$cval, digits = digits)
  }
	cat("\n")
	if(length(x$covariance))
	{
		cat("\n")
		if(length(x$confint))
		{
			cat(gettextf("%d%% Confidence Intervals:",(100*x$level)),"\n")
			print.default(format(x$confint, digits = digits), print.gap = 3,
 			quote = FALSE)
		}
		else cat("No confidence intervals\n")
		cat("\n")
		if(length(x$bootconfint))
		{
			cat(gettextf("%d%% Confidence Intervals (Percentile):", (100*x$level)),"\n")
			print.default(format(x$bootconfint, digits = digits), print.gap = 3,
 			quote = FALSE)
		}
	}
	else
		cat("No covariance matrix has been estimated, hence no t-tests or confidence intervals are returned.\n To get these, choose covar=TRUE in the function call for qme().") 
	cat("\n")
	invisible(x)

}

print.qme <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(x$coefficients))	 
  {
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  cat("Iterations:\n")
  print.default(format(x$counts, digits = digits,widht=2), print.gap = 2,
                quote = FALSE)
  cat("\n")
  cat("Threshold information:\n")
  print.data.frame(x$cval, digits = digits)
  cat("\n")
  invisible(x)
}	


summary.qme <- function(object,level=0.95,...)
{
	x <- object
	if(length(x$covariance))
	{	
		R <- x$covariance
		se <- sqrt(diag(R))
		est <- t(x$coefficients)
		tval <- est/se
		rdf <- x$df.residual
		ans <- list()
		ans$call <- x$call
		ans$counts <- x$counts
    ans$cval <- x$cval
		ans$level <- level   
		ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
		dimnames(ans$coefficients) <- list(rownames(t(x$coefficients)), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
		ans$covariance <- R
		error <- qt((1-(1-level)/2),df=rdf)*se	
		left <- t(x$coefficients)-error
		right <- t(x$coefficients)+error
		ans$confint	<- cbind(left,right)
		dimnames(ans$confint) <- list(rownames(t(x$coefficients)), c("Lower","Upper"))
		ans$bootconfint <- bootconfintQME(x$bootrepl,level)
		dimnames(ans$bootconfint) <- list(rownames(t(x$coefficients)), c("Lower","Upper"))
	}
	else
	{
		est <- t(x$coefficients)
		ans <- list()
		ans$call <- x$call
    ans$counts <- x$counts 
    ans$cval <- x$cval
		se <- NA
		tval <- NA
		pval <- NA  
		ans$coefficients <- cbind(est, se, tval, pval)
		dimnames(ans$coefficients) <- list(rownames(t(x$coefficients)), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
	}
	class(ans) <- c("summary.qme","qme")
	ans

}


coef.qme <- function(object,...)
{
	object$coefficients
}

vcov.qme <- function(object,...) 
{
	if(length(object$covariance))
		object$covariance
	else
		cat("No covariance matrix has been estimated. To do this, choose covar=TRUE in the function call for qme(). \n")

}


residuals.qme <- function(object,...) 
{
	object$residuals
}


fitted.qme <- function(object,...) 
{
	object$fitted.values
}

setMethod("print",signature("qme"),print.qme)

setMethod("summary",signature(object="qme"),summary.qme)

setMethod("coef", signature(object="qme"), coef.qme)

setMethod("vcov", signature(object="qme"), vcov.qme)

setMethod("residuals", signature(object="qme"), residuals.qme)

setMethod("fitted", signature(object="qme"), fitted.qme)

setMethod("print",signature("summary.qme"),print.summary.qme)







