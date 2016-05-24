abc2xy <-
function(abc)
{
	tmat = c(1 / sqrt(3), 0, 2 / sqrt(3), 1, 0, 0)
	dim(tmat) = c(3, 2)
	res = as.matrix(abc) %*% tmat
	colnames(res) = c("x", "y")
	res
}

