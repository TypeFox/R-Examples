"scan.haplo" <- 
		function(formula,data,snpsubset,idsubset,n.slide=2,bcast=10,simulate=FALSE,trait.type,...) 
{
	if (!require(haplo.stats)) {
		stop("this function requires 'haplo.stats' package to be installed")
	}
	if (!is(data,"gwaa.data")) stop("wrong type of data argument, must be gwaa.data")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	
	avrs <- all.vars(as.formula(formula))
	if (!any(avrs=="CRSNP")) stop("formula must contain CRSNP variable to be replaced with the analysis SNPs")
	avrs <- avrs[avrs!="CRSNP"]
#	attach(data@phdata,pos=2,warn.conflicts=FALSE)
	cov <- NA
	if (length(avrs)>1) {
		cov <- phdata(data)[[avrs[2]]]
		if (length(avrs)>2) for (i in 3:length(avrs)) cov <- cbind(cov,phdata(data)[[avrs[i]]])
	}
	tra <- phdata(data)[[avrs[1]]]
	if (missing(trait.type)) {
		if (length(unique(tra))==2) 
			trait.type<-"binomial" 
		else 
			trait.type<-"gaussian"
	}
#	detach(data@phdata)
	
	if (length(unique(tra))<2) stop("Trait is monomorphic!") 
	
	oldmap <- data@gtdata@map
	nsnps <- data@gtdata@nsnps
	Pv <- rep(1,nsnps-n.slide+1)
	name <- rep("",nsnps-n.slide+1)
	map <- rep(-1,nsnps-n.slide+1)
	for (i in 1:(nsnps-n.slide+1)) {
		cumm <- oldmap[i]
		name[i] <- data@gtdata@snpnames[i];
		for (j in (i+1):(i+n.slide-1)) {
			name[i]=paste(name[i],data@gtdata@snpnames[j],sep="-");
			cumm <- cumm + oldmap[j]
		}
		map[i] <- cumm/n.slide
	}
	family <- "haplo.score.slide test"
	for (j in 1:(nsnps-n.slide+1)) {
		hsdta <- as.hsgeno.snp.data(data@gtdata[,j:(j+n.slide-1)])
		tmpo <- haplo.score.slide(y=tra,geno=hsdta,n.slide=n.slide,simulate=simulate,x.adj=cov,trait.type=trait.type,...)
		if (simulate) 
			Pv[j] <- tmpo$df$global.p.sim
		else 
			Pv[j] <- tmpo$df$score.global.p
		if (j/bcast == round(j/bcast)) {
			cat("\b\b\b\b\b\b",round(100*j/(nsnps-n.slide+1),digits=2),"%",sep="")
			flush.console();
		}
	}
	if (bcast<=(nsnps-n.slide+1)) cat("\n")
	#results <- P1df=Pv; snpnames=name,map=map,chromosome=data@gtdata@chromosome[1:(nsnps-n.slide+1)]
	#out <- new("scan.gwaa",results=results,formula=match.call(),family=family,
	#		idnames=idnames(data))
	#class(out) <- "scan.gwaa"
	#out
	results <- data.frame(N=NA,
			effB = NA, se_effB = NA, chi2.1df = NA, P1df = Pv, 
			Pc1df = NA,
			stringsAsFactors = FALSE)
	rownames(results) <- snpnames(data)[1:(nsnps-n.slide+1)]
	out <- new("scan.gwaa",
			results=results,
			annotation = annotation(data[,1:(nsnps-n.slide+1)]), 
			lambda = list(0),
			idnames = idnames(data), 
			call = match.call(), 
			family = family
	) 
	
}
