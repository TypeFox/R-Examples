##
##  f i n d i n t e r v a l s . R  Find Interval Indices
##


findintervals <- function(x, xs) {
    if (length(x) == 0 || length(xs) == 0) return(c())
    if (!is.vector(x, mode="numeric") || !is.vector(x, mode="numeric"))
        stop("Arguments 'x' and 'xs' must be numeric vectors.")
    if (length(x) != 1)
        stop("Length of vector 'x' must be 1.")

    n <- length(xs)
    xsx <- xs - x

	i0 <- which(xsx == 0)
	i1 <- which(xsx[1:(n-1)] * xsx[2:n] < 0)
	
	return(sort(c(i0, i1)))
}
