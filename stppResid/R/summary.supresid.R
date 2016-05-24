summary.supresid <- function(object, ...)
{
  x <- object
	k <- x$k
	n <- nrow(x$residuals)
	vol <- diff(x[[1]]$xcoord) * diff(x[[1]]$ycoord) * diff(x[[1]]$tcoord)
	n.exp <- k * vol
	if(n < n.exp)
		p.val <- ppois(n, n.exp) 
	if(n > n.exp)
		p.val <- ppois(n, n.exp, lower.tail = FALSE) 
	if(n == n.exp)
		p.val <- 1
	Y <- list(k = k, n = n, n.exp = n.exp, p.val = p.val)
	class(Y) <- "summary.supresid"
	return(Y)
}