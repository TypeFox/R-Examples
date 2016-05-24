duplicated.matrix = function (x, incomparables = FALSE, MARGIN = 1L, fromLast = FALSE, ...)
{
	if (!is.matrix(x) ||!is.numeric(x[1]) || !identical(incomparables, FALSE) || MARGIN!=1L)
		return(base::duplicated.matrix(x, incomparables, MARGIN, fromLast, ...))
	.Call(C_dupRowNumMat,x, as.logical(fromLast))
}

unique.matrix=function (x, incomparables = FALSE, MARGIN = 1, fromLast = FALSE, ...)
{
	if (!is.matrix(x) ||!is.numeric(x[1]) || !identical(incomparables, FALSE) || MARGIN!=1L)
		return(base::unique.matrix(x, incomparables, MARGIN, fromLast, ...))
	x[!.Call(C_dupRowNumMat,x,as.logical(fromLast)),,drop=FALSE]
}
