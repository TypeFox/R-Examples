#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

confint.gel <- function(object, parm, level = 0.95, lambda = FALSE, ...)
		{
		z <- object	
		n <- nrow(z$gt)
		
		se_par <- sqrt(diag(z$vcov_par))
		par <- z$coefficients
		tval <- par/se_par

		se_parl <- sqrt(diag(z$vcov_lambda))
		lamb <- z$lambda

		zs <- qnorm((1 - level)/2, lower.tail=FALSE)
		ch <- zs*se_par

		if(!lambda)
			{
			ans <- cbind(par-ch, par+ch)
			dimnames(ans) <- list(names(par), c((1 - level)/2, 0.5+level/2))
			}
		if(lambda)
			{
			chl <- zs*se_parl
			ans <- cbind(lamb - chl, lamb + chl)
			dimnames(ans) <- list(names(lamb), c((1 - level)/2, 0.5 + level/2))
			}		
		if(!missing(parm))
			ans <- ans[parm,]
		ans
		}

coef.gel <- function(object, lambda = FALSE, ...) 
	{
	if(!lambda)
		object$coefficients
	else
		object$lambda
	}

vcov.gel <- function(object, lambda = FALSE, ...) 
	{
	if(!lambda)
		object$vcov_par
	else
		object$vcov_lambda
	}

print.gel <- function(x, digits = 5, ...)
	{
	if (is.null(x$CGEL))
		cat("Type de GEL: ", x$type, "\n")
	else
		cat("CGEL of type: ", x$type, " (alpha = ", x$CGEL, ")\n")
	if (!is.null(attr(x$dat,"smooth")))
		{
		cat("Kernel: ", attr(x$dat,"smooth")$kernel," (bw=",
		attr(x$dat,"smooth")$bw,")\n\n")
		}
	else
		cat("\n")

	cat("Coefficients:\n")
	print.default(format(coef(x), digits = digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	cat("Lambdas:\n")
	print.default(format(coef(x, lambda = TRUE), digits = digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	cat("Convergence code for the coefficients: ", x$conv_par,"\n")
	cat("Convergence code for Lambda: ", x$conv_lambda$convergence,"\n")
	invisible(x)
	}

print.summary.gel <- function(x, digits = 5, ...)
	{
	cat("\nCall:\n")
	cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	if (is.null(x$CGEL))
		cat("Type of GEL: ", x$type, "\n")
	else
		cat("CGEL of type: ", x$type, " (alpha = ", x$CGEL, ")\n")

	if (!is.null(x$smooth))
		{
		cat("Kernel: ", x$smooth$kernel," (bw=", x$smooth$bw,")\n\n")
		}
	else
		cat("\n")

	cat("Coefficients:\n")
	print.default(format(x$coefficients, digits = digits),
                      print.gap = 2, quote = FALSE)

	cat("\nLambdas:\n")
	print.default(format(x$lambda, digits=digits),
                      print.gap = 2, quote = FALSE)

	cat("\n", x$stest$ntest, "\n")
	print.default(format(x$stest$test, digits=digits),
                      print.gap = 2, quote = FALSE)

	cat("\nConvergence code for the coefficients: ", x$conv_par, "\n")
	cat("\nConvergence code for the lambdas: ", x$conv_lambda$convergence, "\n")
	
	invisible(x)
	}

summary.gel <- function(object, ...)
	{
	z <- object
	n <- nrow(z$gt)
	se_par <- sqrt(diag(z$vcov_par))
	par <- z$coefficients
	tval <- par/se_par

	se_parl <- sqrt(diag(z$vcov_lambda))
	lamb <- z$lambda
	tvall <- lamb/se_parl

	ans <- list(type = z$type, call = z$call)
	names(ans$type) <-"Type of GEL"
	
	ans$coefficients <- round(cbind(par, se_par, tval, 2 * pnorm(abs(tval), lower.tail = FALSE)), 5)
	ans$lambda <- round(cbind(lamb,se_parl, tvall, 2 * pnorm(abs(tvall), lower.tail = FALSE)), 5)

    	dimnames(ans$coefficients) <- list(names(z$coefficients), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    	dimnames(ans$lambda) <- list(names(z$lambda), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

	ans$stest=specTest(z)

	if (z$type == "EL")
		ans$badrho <- z$badrho
	if (!is.null(z$weights))
		{
		ans$weights <- z$weights
		}
	ans$conv_par <- z$conv_par
	ans$conv_pt <- z$conv_pt
	ans$conv_moment <- cbind(z$conv_moment)
	ans$conv_lambda <- z$conv_lambda
	ans$CGEL <- z$CGEL
	if (!is.null(attr(object$dat,"smooth")))
		ans$smooth <- attr(object$dat,"smooth")
	names(ans$conv_pt) <- "Sum_of_pt"
	dimnames(ans$conv_moment) <- list(names(z$gt), "Sample_moment_with_pt")
	class(ans) <- "summary.gel"
	ans	
}

residuals.gel <- function(object, ...) 
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$residuals
	}

fitted.gel <- function(object, ...)
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$fitted.value
	}

formula.gel <- function(x, ...)
{
    if(is.null(x$terms))
	stop("The gel object was not created by a formula")
    else
	formula(x$terms)
}

estfun.gel <- function(x, ...)
  {
  stop("estfun is not yet available for gel objects")
  }

bread.gel <- function(x, ...)
  {
  stop("Bread is not yet available for gel objects")
  }



