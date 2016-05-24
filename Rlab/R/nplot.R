"nplot" <-

function(x, ylab = "Quantiles of Standard Normal", xlab = deparse(



	substitute(x)), full = TRUE, ...)



{



	out <- is.na(x)



	if(any(out)) {



		x <- x[!out]



	}



	n <- length(x)



	if(full) {



		y <- qnorm((rank(x) - 0.5)/n)



	}



	else {



		x <- abs(x)



		y <- qnorm((rank(x) - 0.5)/n)



	}



	plot(x, y, xlab = xlab, ylab = ylab, ...)



	if(is.null(names(x)))



		names(x) <- format(1.:length(x))



	text(x, y, names(x), adj = 0.)



}

