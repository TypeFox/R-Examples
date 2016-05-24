msBP.nrvTrees <- function(sh, maxS = max(sh[,1]))
{
	N <- nrow(sh)
	veclen <- 2^(maxS+1) - 1
	empty <- rep(0, veclen)
	res <- .C("allTrees_C", as.integer(sh[,1]), as.integer(sh[,2]), as.integer(maxS), as.integer(N), n=as.double(empty), r=as.double(empty), v=as.double(empty), PACKAGE = "msBP")
	n <- vec2tree(res$n)
	r <- vec2tree(res$r)
	v <- vec2tree(res$v)
	list(n=n, r=r, v=v)
}


