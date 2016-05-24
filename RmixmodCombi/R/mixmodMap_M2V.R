mixmodMap_M2V <-
function(M, n = nrow(M), K = ncol(M))
{
	z = apply(M, 1, function(x) which(as.logical(x)))
	z
}
