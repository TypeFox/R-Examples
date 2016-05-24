xy2abc <-
function(xy)
{
	tmat = c(sqrt(3) / 2, -0.5, 0, 1)
	dim(tmat) = rep(2, 2)
	res = as.matrix(xy) %*% tmat
	res = cbind(res[,2], 1 - rowSums(res), res[,1])
	colnames(res) = c("a", "b", "c")
	res
}

