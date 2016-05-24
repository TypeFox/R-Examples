ExoLVsCor <-
function(N,n,bex,rie,y.ex,a.nle,a.ie)
	{
	y.ex.tot <- array(matrix(0,n,bex+prod(c(1:ncol(rie)))),dim=c(n,bex+prod(c(1:ncol(rie))),N))
	y.ex.tot <- abind(y.ex,a.nle,a.ie,along=2)

	y.ex.cor <- array(matrix(0,bex+prod(c(1:ncol(rie))),bex+prod(c(1:ncol(rie)))),dim=c(bex+prod(c(1:ncol(rie))),bex+prod(c(1:ncol(rie))),N))
	for (i in 1:N)
		{
		y.ex.cor[,,i] <- cov(y.ex.tot[,,i]) 
		}
	list(y.ex.tot=y.ex.tot,y.ex.cor=y.ex.cor)
}

