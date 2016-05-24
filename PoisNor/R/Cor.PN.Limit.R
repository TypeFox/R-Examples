Cor.PN.Limit <-
function(lam){
	X=rnorm(100000,0,1)
	Y=rnorm(100000,0,1)

	U = pnorm(X)

	Xpois = qpois(U,lam)

	max = cor(Xpois[order(Xpois)],Y[order(Y)])
	min = cor(Xpois[order(Xpois,decreasing=TRUE)],Y[order(Y)])
	rm(X,Y)
	return(c(min,max))
}
