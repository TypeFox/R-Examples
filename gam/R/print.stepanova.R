"print.stepanova"<-
function(x, digits = .Options$digits, quote = F, drop = F, ...)
{
	heading <- attr(x, "heading")
	if(!is.null(heading))
		cat(heading, sep = "\n")
	attr(x, "heading") <- NULL
	d <- dim(x)
	for(i in 1:d[2]) {
		xx <- x[[i]]
		if(!length(levels(xx)) && is.numeric(xx)) {
			xna <- is.na(xx)
			xx <- format(zapsmall(xx, digits))
			xx[xna] <- ""
			x[[i]] <- xx
		}
	}
	if(d[1] == 1 && drop) {
		x <- t(as.matrix(x))
		dn <- dimnames(x)
		dn <- paste(dn[[1]], ":", sep = "")
		dimnames(x) <- list(dn, "")
	}
	NextMethod("print")
}
