dev2div <- function(x, n)  {
	if (class(x)[[1]]=="mdm") exp(x$deviance/2/nrow(x$fitted.values))
	else exp(x/2/n)
}

