hdp <-
function(x,y)
{
	x <- as.matrix(x)
	y <- as.matrix(y)
	if(nrow(x) != nrow(y))stop("row dimensions are not same")
	z <- x[,rep(1:ncol(x),rep(ncol(y),ncol(x)))]*y[,rep(1:ncol(y),ncol(x))]
	z
}
