"kurtosis" <-
function (x, na.rm = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, kurtosis, na.rm = na.rm)
    else if (is.vector(x)) {
	if (na.rm) x <- x[!is.na(x)] 
	n <- length(x)
	n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2)
	}
    else if (is.data.frame(x)) 
        sapply(x, kurtosis, na.rm = na.rm)
    else kurtosis(as.vector(x), na.rm = na.rm)
}

