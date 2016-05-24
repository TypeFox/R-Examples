"mmscore" <-
		function(h2object,data,snpsubset,idsubset,strata,times=1,quiet=FALSE,
				bcast=10,clambda=TRUE,propPs=1.0) {
	if (is(data,"gwaa.data")) 
	{
		checkphengen(data)
		data <- data@gtdata
	}
	if (!is(data,"snp.data")) {
		stop("wrong data class: should be gwaa.data or snp.data")
	}
	if (class(h2object) != "polygenic") stop("h2object must be of polygenic-class")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (missing(strata)) {nstra=1; strata <- rep(0,data@nids)}
	
	if (length(strata)!=data@nids) stop("Strata variable and the data do not match in length")
	if (any(is.na(strata))) stop("Strata variable contains NAs")
		
	tmeas <- h2object$measuredIDs
	resid <- h2object$residualY
	if (any(tmeas == FALSE)) {
		if (!quiet) warning(paste(sum(!tmeas),"people (out of",length(tmeas),") excluded because they have trait or covariate missing\n"),immediate. = TRUE)
		if (length(tmeas) != data@nids) stop("Dimension of the outcome and SNP data object are different")
		data <- data[tmeas,]
		strata <- strata[tmeas]
		resid <- resid[tmeas]
	}
	if (any(strata!=0)) {
		olev <- levels(as.factor(strata))
		nstra <- length(olev)
		tstr <- strata
		for (i in 0:(nstra-1)) tstr <- replace(tstr,(strata==olev[i+1]),i)
		strata <- tstr
		rm(tstr)
	}
	nstra <- length(levels(as.factor(strata)))
	
	lenn <- data@nsnps;
	out <- list()
	for (j in c(1:(times+1*(times>1)))) {
		if (j>1) resid <- sample(resid,replace=FALSE)
		chi2 <- .C("mmscore_20110916",as.raw(data@gtps),as.double(resid),as.double(h2object$InvSigma),as.integer(data@nids),as.integer(data@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(7*data@nsnps), PACKAGE="GenABEL")$chi2
		if (any(data@chromosome=="X")) {
			ogX <- data[,data@chromosome=="X"]
			sxstra <- strata; sxstra[ogX@male==1] <- strata[ogX@male==1]+nstra
			chi2X <- .C("mmscore_20110916",as.raw(ogX@gtps),as.double(resid),as.double(h2object$InvSigma),as.integer(ogX@nids),as.integer(ogX@nsnps), as.integer(nstra*2), as.integer(sxstra), chi2 = double(7*ogX@nsnps), PACKAGE="GenABEL")$chi2
			revec <- (data@chromosome=="X")
			revec <- rep(revec,6)
			chi2 <- replace(chi2,revec,chi2X)
			rm(ogX,chi2X,revec);gc(verbose=FALSE)
		}
		if (j == 1) {
			chi2.1df <- chi2[1:lenn];
			chi2.2df <- chi2[(lenn+1):(2*lenn)];
			#out$chi2.1df <- chi2.1df
			#out$chi2.2df <- chi2.2df
			actdf <- chi2[(2*lenn+1):(3*lenn)];
			lambda <- list()
			if (is.logical(clambda)) {
				if (lenn<10) {
					warning("no. observations < 10; Lambda set to 1")
					lambda$estimate <- 1.0
					lambda$se <- NA
				} else {
					if (lenn<100) warning("Number of observations < 100, Lambda estimate is unreliable")
					lambda <- estlambda(chi2.1df,plot=FALSE,proportion=propPs)
					if (lambda$estimate<1.0 && clambda==TRUE) {
						warning("Lambda estimated < 1, set to 1")
						lambda$estimate <- 1.0
						lambda$se <- NA
					}
				}
			} else {
				if (is.numeric(clambda)) {
					lambda$estimate <- clambda
					lambda$se <- NA
				} else if (is.list(clambda)) {
					if (any(is.na(match(c("estimate","se"),names(clambda)))))
						stop("when clambda is list, should contain estimate and se")
					lambda <- clambda
					lambda$se <- NA
				} else {
					stop("clambda should be logical, numeric, or list")
				}
			}
			chi2.c1df <- chi2.1df/lambda$estimate
			effB <- chi2[(3*lenn+1):(lenn*4)]
#			effAB <- chi2[(4*lenn+1):(lenn*5)]
#			effBB <- chi2[(5*lenn+1):(lenn*6)]
			if (times>1) {
				pr.1df <- rep(0,lenn)
#				pr.2df <- rep(0,lenn)
				pr.c1df <- rep(0,lenn)
			}
		} else {
			th1 <- max(chi2[1:lenn])
			pr.1df <- pr.1df + 1*(chi2.1df < th1)
#			pr.2df <- pr.2df + 1*(chi2.2df < max(chi2[(lenn+1):(2*lenn)]))
			pr.c1df <- pr.c1df + 1*(chi2.c1df < th1)
			if (!quiet && ((j-1)/bcast == round((j-1)/bcast))) {
				cat("\b\b\b\b\b\b",round((100*(j-1)/times),digits=2),"%",sep="")
				flush.console()
			}
		}
	}
	if (times > bcast) cat("\n")
	
	if (times>1) {
		P1df <- pr.1df/times
		P1df <- replace(P1df,(P1df==0),1/(1+times))
#		out$P2df <- pr.2df/times
#		out$P2df <- replace(out$P2df,(out$P2df==0),1/(1+times))
		Pc1df <- pr.c1df/times
#		out$Pc1df <- replace(out$Pc1df,(out$Pc1df==0),1/(1+times))
	} else {
		P1df <- pchisq(chi2.1df,1,lower.tail=F)
#		out$P2df <- pchisq(chi2.2df,actdf,lower.tail=F)
		Pc1df <- pchisq(chi2.c1df,1,lower.tail=F)
	}
	#out$lambda <- lambda
	#out$effB <- effB #*var(h2object$residualY,na.rm=T)
#	out$effAB <- effAB
#	out$effBB <- effBB
	#out$snpnames <- data@snpnames
	#out$map <- data@map
	#out$chromosome <- data@chromosome
	#out$idnames <- data@idnames
	#out$formula <- match.call()
	#out$family <- paste("score test for association with trait type") #,trait.type)
	effAB <- rep(NA,length(P1df))
	effBB <- rep(NA,length(P1df))
	P2df <- rep(NA,length(P1df))
	#out$N <- chi2[(6*lenn+1):(lenn*7)]
	#class(out) <- "scan.gwaa"
	#out
	if (is.null(Pc1df)) {
		results <- data.frame(N=chi2[(6*lenn+1):(lenn*7)],
				effB = effB, se_effB = abs(effB)/sqrt(chi2.1df), chi2.1df = chi2.1df, P1df = P1df, 
				effAB=effAB, effBB=effBB, chi2.2df = chi2.2df, P2df = P2df,
				stringsAsFactors = FALSE)
	} else {
		results <- data.frame(N=chi2[(6*lenn+1):(lenn*7)],
				effB = effB, se_effB = abs(effB)/sqrt(chi2.1df), chi2.1df = chi2.1df, P1df = P1df, 
				Pc1df = Pc1df, 
				effAB=effAB, effBB=effBB, chi2.2df = chi2.2df, P2df = P2df,
				stringsAsFactors = FALSE)
	}
	rownames(results) <- snpnames(data)
	out <- new("scan.gwaa",
			results=results,
			annotation = annotation(data), 
			lambda = lambda,
			idnames = idnames(data), 
			call = match.call(), 
			family = "mmscore"
	) 
	out
}
