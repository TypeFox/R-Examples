"formetascore_engage" <-
function(formula,data,stat=qtscore,transform=ztransform, build="unknown", ...) {
	if (is.character(transform)) {
		if (transform=="no") trvar <- formula
		else stop("transform argument not recognised")
	} else {
		trvar <- transform(formula=formula,data=data)
	}
	data <- data[!is.na(trvar),]
	trvar <- trvar[!is.na(trvar)]
	res <- stat(trvar,data,...)
	sum <- summary(data@gtdata)
	callr <- sum$Call
	Phwe <- sum$Pexact
	efff <- sum$Q.2
	reff <- 1. - sum$Q.2
#	rm(sum);gc(verbose=FALSE)
	serr <- res$effB/sqrt(res$chi2.1df)
	coding <- as.character(data@gtdata@coding)
 	out <- data.frame(name=res$snpnames,chromosome=res$chromosome,
		position=res$map,strand=as.character(data@gtdata@strand),
#		coding=coding,
		allele1=alleleID.reference()[coding],
		allele2=alleleID.effective()[coding],
		build=rep(build,length(coding)),
		effallele=alleleID.effective()[coding],
		effallelefreq=efff,
		n=res$N,beta=res$effB,sebeta=serr,p=res$P1df,
		pgc=res$Pc1df,lambda=rep(res$lam$est,length(coding)),
		pexhwe=Phwe,call=callr,fmax=sum$Fmax,plrthwe=sum$Plrt,stringsAsFactors=F)
	rownames(out) <- res$snpnames
	out
}

