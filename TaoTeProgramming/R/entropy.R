"entropy" <-
function (xmat) 
{
	ent <- function(x) length(rle(x)$lengths)
	entrow <- apply(xmat, 1, ent)
	entcol <- apply(xmat, 2, ent)
	sum(c(entrow, entcol))
}

