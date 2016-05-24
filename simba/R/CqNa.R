"CqNa" <-
function(x, q = 0, base = exp(1), ...){
#	if(q==1){ stop("The function is currently not defined at q = 1") }
	mat <- x
	n <- nrow(mat)
	e <- q-1
	# sampling unit weights (needed for relative homogeneity)
	su.w <- rowSums(mat/sum(mat))
	# true diversity of the weights (needed for relative homogeneity)
	td.w <- exp(-sum(su.w*log(su.w, base=base), na.rm=TRUE))
	# true beta via trudi()
	td.b <- trudi(mat, q=q, ...)[2]
	if(q!=1){
		# cqn (overlap of order q, not defined for q=1)
		cqn <- ((1/td.b)^e - (1/n)^e) / (1 - (1/n)^e)
	}
	if(q==1){
		hb.shan <- as.numeric(diff(log(trudi(mat, q=1))[c(3,1)]))
		cqn <- (log(n, base = base) - hb.shan) / log(n, base = base)
	}
	# homogeniety sensu MacArthur but ranged between 0 and 1
	hom.MA <- (1/td.b - 1/n) / (1 - 1/n)
	# relative homogeneity (homogeneity sensu MacArthur but
	# for q=1 and ranged between 0 an 1 by using the true diversity of the weights)
	hom.rel <- (1/td.b - 1/td.w) / (1 - 1/td.w)
	res <- c(cqn=cqn, hom.MA=hom.MA, hom.rel=hom.rel)
	names(res)[1] <- paste("C", q, n, "o", sep="")
	return(res)
}