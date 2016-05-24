"formetascore" <-
		function(formula,data,stat=qtscore,transform="no",build="unknown",verbosity=1, ...) {
	if (!is(data,"gwaa.data")) stop("data argument must have gwaa.data-class")
	checkphengen(data)
	
	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		formula <- phdata(data)[[as(match.call()[["formula"]],"character")]] 
	}
	
	#	if (!missing(data)) attach(data@phdata,pos=2,warn.conflicts=FALSE)
	if (is(formula,"polygenic")) {
		pm <- pmatch("stat",names(match.call()))
		pm <- (pm[!is.na(pm)])[1]
		a <- match.call()[[pm]]
		if (a != "mmscore" && a != "mmscore()") stop("stat should be mmscore when polygenic object is analysed")
		if (!is(transform,"character")) stop("transform should be \"no\" when polygenic object is analysed")
		if (transform != "no") stop("transform should be \"no\" when polygenic object is analysed")
	}
	if (is.character(transform)) {
		if (transform!="no") stop("transform argument not recognised")
	} else {
		formula <- transform(formula=formula,data=data)
	}
	if (is(formula,"formula")) {
		mf <- model.frame(formula,data@phdata,na.action=na.omit,drop.unused.levels=TRUE)
		mids <- rownames(data@phdata) %in% rownames(mf)
	} else if (is(formula,"polygenic")) {
		mids <- which(!is.na(formula$residualY))
	} else {
		mids <- which(!is.na(formula))
	}
#	if (!missing(data)) detach(data@phdata)
	if (verbosity<0) stop("verbosity parameter must be positive integer")
	res <- stat(formula,data,...)
	sum <- summary(data@gtdata[mids,])
	callr <- sum$Call
	Phwe <- sum$Pexact
	efff <- sum$Q.2
	reff <- 1. - sum$Q.2
	serr <- abs(res[,"effB"])/sqrt(res[,"chi2.1df"])
	#coding <- coding(data)
	if (verbosity==0)
		out <- data.frame(
				name=snpnames(res),chromosome=chromosome(res),
				position=map(res),strand=strand(res),
				allele1=refallele(res),
				allele2=effallele(res),
				effallele=effallele(res),
				beta=res[,"effB"],sebeta=res[,"se_effB"],p=res[,"P1df"],
				stringsAsFactors=FALSE
		)
	else if (verbosity==1)
		out <- data.frame(
				name=snpnames(res),chromosome=chromosome(res),
				position=map(res),strand=strand(res),
				allele1=refallele(res),
				allele2=effallele(res),
				effallele=effallele(res),
				effallelefreq=efff,
				n=res[,"N"],
				beta=res[,"effB"],sebeta=res[,"se_effB"],p=res[,"P1df"],
				pgc=res[,"Pc1df"],
				pexhwe=Phwe,call=callr,stringsAsFactors=FALSE
		)
	else 
		out <- data.frame(
				name=snpnames(res),chromosome=chromosome(res),
				position=map(res),strand=strand(res),
				allele1=refallele(res),
				allele2=effallele(res),
				build=rep(build,length(coding)),
				effallele=effallele(res),
				effallelefreq=efff,
				n=res[,"N"],
				beta=res[,"effB"],sebeta=res[,"se_effB"],p=res[,"P1df"],
				pgc=res[,"Pc1df"],lambda=lambda(res),
				pexhwe=Phwe,call=callr,stringsAsFactors=FALSE
		)
	
	
	rownames(out) <- snpnames(res)
	out
}
