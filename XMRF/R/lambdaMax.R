lambdaMax <-
function(X){
	tmp = t(X)%*%X
	return(1/nrow(X)*max( abs(tmp[upper.tri(tmp)])))
}
