"scan.haplo.2D" <- 
function(formula,data,snpsubset,idsubset,bcast=10,simulate=FALSE,trait.type,...) {

	if (!require(haplo.stats)) {
		stop("this function requires 'haplo.stats' package to be installed")
	}
	if (!is(data,"gwaa.data")) stop("wrong type of data argument, must be gwaa.data")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]

	avrs <- all.vars(as.formula(formula))
	if (!any(avrs=="CRSNP")) stop("formula must contain CRSNP variable to be replaced with the analysis SNPs")
	avrs <- avrs[avrs!="CRSNP"]
	#attach(data@phdata,pos=2,warn.conflicts=FALSE)
	cov <- NA
	if (length(avrs)>1) {
		cov <- matrix(phdata(data)[[avrs[2]]],ncol=1)
		if (length(avrs>2)) for (i in 3:length(avrs)) cov <- cbind(cov,phdata(data)[[avrs[i]]])
	}
	tra <- phdata(data)[[avrs[1]]]
	if (missing(trait.type)) {
		if (length(unique(tra))==2) 
			trait.type<-"binomial" 
		else 
			trait.type<-"gaussian"
	}
	#detach(data@phdata)

        if (length(unique(tra))<2) stop("Trait is monomorphic!") 

	nsnps <- data@gtdata@nsnps
        Pv <- matrix(rep(NA,(nsnps*nsnps)),nrow=nsnps)
	formula <- match.call()
	family <- "haplo.score.slide test"
	donan <- 0
	for (j1 in 1:(nsnps-1)) {
	for (j2 in (j1+1):nsnps) {
	  twonam <- c(data@gtdata@snpnames[j1],data@gtdata@snpnames[j2])
	  hsdta <- as.hsgeno.snp.data(data@gtdata[,twonam])
	  tmpo <- haplo.score.slide(y=tra,geno=hsdta,n.slide=2,simulate=simulate,trait.type=trait.type,x.adj=cov,...)
	if (simulate) {
		Pv[j1,j2] <- tmpo$df$global.p.sim
		if (Pv[j1,j2] <= 0) Pv[j1,j2] <- 1./(score.sim.control()$max.sim+1.)
	}
	else 
	  	Pv[j1,j2] <- tmpo$df$score.global.p
	  donan <- donan + 1
	  if (bcast) {
	    if (donan/bcast == round(donan/bcast)) {
		  cat("\b\b\b\b\b\b",round(100*(donan/((nsnps*(nsnps-1))/2)),digits=2),"%",sep="")
		  flush.console()
	    }
	  }
	}
	}
        if (bcast<100.01) cat("\n")
	rownames(Pv) <- data@gtdata@snpnames
	colnames(Pv) <- data@gtdata@snpnames 
  	Pint1df <- rep(NA,length(Pv))
  	Pint2df <- Pint1df
	out <- list(P1df=Pv,Pint1df=Pint1df,Pint2df=Pint2df,P2df=Pv,snpnames=data@gtdata@snpnames,formula=match.call(),family=family,map=data@gtdata@map,idnames=data@gtdata@idnames)
  	out$Pc1df <- rep(NA,length(Pv))
	class(out) <- "scan.gwaa.2D"
	out
}
