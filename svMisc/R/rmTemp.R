rmTemp <- function (x)
{
	if (!is.character(x))
		stop("'x' must be character string(s)!")
	l <- length(x)
	res <- rep(TRUE, l)
	if (l > 1) names(res) <- x
	for (i in 1:l) {
		exst <- exists(x[i], envir = TempEnv(), inherits = FALSE)
		res0 <- try(if (exst) rm(list = x[i],
			envir = TempEnv()), silent = TRUE)
		if (!exst || inherits(res0, "try-error")) res[i] <- FALSE
	}
	invisible(res)
}
