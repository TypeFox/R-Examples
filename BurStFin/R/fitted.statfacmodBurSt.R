"fitted.statfacmodBurSt" <-
function (object, output="full", ...)
{
        fun.copyright <- "Placed in the public domain 2006-2009 by Burns Statistics"
	fun.version <- "fitted.statfacmodBurSt 005"

	if(!is.character(output) || length(output) != 1) {
		stop(paste("'output' should be a single character string",
			"-- given has mode", mode(output), "and length",
			length(output)))
	}
	output.menu <- c("full", "systematic", "specific")
	output.num <- pmatch(output, output.menu, nomatch=0)
	if(output.num == 0) {
		stop(paste("unknown or ambiguous input for 'output'",
			"-- the allowed choices are:",
			paste(output.menu, collapse=", ")))
	}
	output <- output.menu[output.num]
	switch(output,
		full={
			ans <- object$loadings %*% t(object$loadings)
			ans <- t(object$sdev * ans) * object$sdev
			diag(ans) <- diag(ans) + object$uniquenesses * 
				object$sdev^2
		},
		systematic={
			ans <- object$loadings %*% t(object$loadings)
			ans <- t(object$sdev * ans) * object$sdev
		},
		specific={
			ans <- diag(object$uniquenesses * object$sdev^2)
			dimnames(ans) <- list(names(object$sdev),
				names(object$sdev))
		}
	)
	attr(ans, "number.of.factors") <- ncol(object$loadings)
	attr(ans, "timestamp") <- object$timestamp
	ans
}

