`invtnp2` <-
function (X, coeff, lengths, lengthsremove, pointsin, lca,nadd = nrow(lca)) 
{

if(is.list(X)){
	coeff<-X$coeff
	lengths<-X$lengths
	lengthsremove<-X$lengthsremove
	lca<-X$lca
	pointsin<-X$pointsin
	X<-X$x
}

if(nadd>0){
	N <- length(pointsin)
	lr <- nrow(lca)
	nc<- ncol(lca)

	nr<-lr-nadd

	ans<-.C("invtnp",as.double(X),coeff=as.double(coeff),lengths=as.double(lengths),
		lengthsremove=as.double(lengthsremove),as.integer(pointsin),
		as.double(t(lca)),as.integer(nadd),as.integer(N),as.integer(lr),
		as.integer(nc),outpo=as.integer(rep(0,times=N+nadd)),
		outlen=as.double(rep(0,times=N+nadd)),PACKAGE="adlift")

	coeff<-ans$coeff
	lengths<-ans$outlen
	pointsin<-ans$outpo
}

return(list(coeff = coeff, lengths = lengths, pointsin = pointsin))
}

