`PointsUpdate2` <-
function (X, coeff, nbrs, index, remove, pointsin, weights, lengths) 
{

r <- which(pointsin == remove)
q <- 0
N <- length(pointsin)
alpha <- matrix(0, 1, length(nbrs))
if (length(nbrs) == 1) {
	q <- which(pointsin == nbrs)
}

ans <- .C("pointsupdate", as.double(X), c = as.double(coeff), 
       	as.integer(length(nbrs)), as.integer(index), as.integer(remove), 
       	as.integer(pointsin), w = as.double(weights), l = as.double(lengths), 
       	N = as.integer(N), a = as.double(alpha), r = as.integer(0), PACKAGE = "adlift")

return(list(coeff = ans$c, lengths = ans$l, r = ans$r, N = ans$N, weights = ans$w, alpha = ans$a))
}

