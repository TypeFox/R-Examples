normalize <-
function(x)
{
	x[is.na(x)] <- 0
	if(is.matrix(x)==T) norm.x <- sweep(x, 2, apply(x, 2, function(k) max(k, na.rm=T)), "/")
	if(is.matrix(x)==F) norm.x <- x/max(x, na.rm=T)
	norm.x[is.na(norm.x)] <- 0
	norm.x
}
