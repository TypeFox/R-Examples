hcr<-function(d, nout, nsampl=1000){
	dmat = as.matrix(d)
	nplots = nrow(dmat)
	if(nout > nplots) stop("nout cannot be larger than the number of plots")
	meand = numeric(nsampl)
	vard = numeric(nsampl)
	sel<-matrix(FALSE,nrow=nout, ncol=nsampl)
	for(s in 1:nsampl) {
		sel[,s] <- sample(x=nplots,size=nout)
		dvec = as.vector(as.dist(dmat[sel[,s],sel[,s]]))
		meand[s] = mean(dvec)
		vard[s] = var(dvec)
	}
	#print(head(sort(meand, decreasing=TRUE)))
	rankdecmean = rank(-meand)
	#print(head(sort(vard)))	
	rankincvar = rank(vard)
	rankfinal = rank(rankdecmean+rankincvar)
	#print(cbind(meand, rankdecmean, vard, rankincvar,rankfinal))
	finalsel = sort(sel[,which.min(rankfinal)])
	return(finalsel)
}