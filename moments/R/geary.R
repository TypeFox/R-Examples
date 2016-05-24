"geary" <-
function (x, na.rm = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, geary, na.rm = na.rm)
    else if (is.vector(x)) {
	if (na.rm) x <- x[!is.na(x)] 
	n <- length(x)
	rho <- sqrt(sum((x-mean(x))^2)/n);
	tau <- sum(abs(x-mean(x)))/n;
	tau/rho
	}
    else if (is.data.frame(x)) 
        sapply(x, geary, na.rm = na.rm)
    else geary(as.vector(x), na.rm = na.rm)
}

