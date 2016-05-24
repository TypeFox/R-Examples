"check.marker.internal" <-
function(data, snpsubset, idsubset,
			callrate=0.95,perid.call=0.95,
			het.fdr=0.01, ibs.threshold = 0.95, ibs.mrk = 2000, ibs.exclude="lower",
			maf, p.level=-1, 
			fdrate = 0.2, hweidsubset, redundant="no", minconcordance = 2.0, 
			qoption="bh95", imphetasmissing = TRUE) {

# qoption = "bh95" (Benjamini & Hochberg 1995) or "storey" (Storey 2003, requires qvalue library)
	if (is(data,"gwaa.data")) {
		if (!missing(snpsubset)) data <- data@gtdata[,snpsubset]
		if (!missing(idsubset)) data <- data@gtdata[idsubset,]
		if (missing(idsubset) & missing(snpsubset)) data <- data@gtdata
	} else if (is(data,"snp.data")) {
		if (!missing(snpsubset)) data <- data[,snpsubset]
		if (!missing(idsubset)) data <- data[idsubset,]
	} else {
		stop("data argument should be of type gwaa.data or snp.data");
	}
	if (missing(maf)) maf <- 5/(2*data@nids)

	nts <- data@nsnps
	cat(nts,"markers and",data@nids,"people in total\n")
	flush.console()
	out <- list()
# check redundancy
#	browser()
	redok <- rep(1,data@nsnps)
	if (redundant != "no" && (redundant == "bychrom" || redundant=="all")) {
		out$details.redundancy <- redundant(data,pairs=redundant, minconcordance = minconcordance)
		out$redundant <- out$details.redundancy[["all"]]
		redok[match(out$details.redundancy$all,data@snpnames)] <- 0
		ntmp <- length(out$redundant); ptmp <- 100*ntmp/nts
		cat(ntmp," (",ptmp,"%) markers excluded as redundant (option = \"",redundant,"\")\n",sep="")
		flush.console()
	} else {
		out$details.redundancy[["all"]] <- NULL;
		out$redundant <- NULL;
	}
# run summary
	s <- summary(data)
	if (any(as.character(data@chromosome) == "Y") && any(data@male==1)) {
		sY <- summary(data[which(data@male==1),(which(as.character(data@chromosome) == "Y"))])
		s[(which(as.character(data@chromosome) == "Y")),] <- sY
	}
# check frequency
	freqok <- pmin(s$Q.2,1-s$Q.2)
	freqok <- (freqok>=maf)
	out$nofreq <- data@snpnames[which(freqok == FALSE)]
	ntmp <- length(out$nofreq); ptmp <- 100*ntmp/nts
	cat(ntmp," (",ptmp,"%) markers excluded as having low (<",maf*100,"%) minor allele frequency\n",sep="")
	flush.console()
# check call rate
	callok <- (s[,"CallRate"]>=callrate)
	out$nocall <- data@snpnames[which(callok == FALSE)]
	ntmp <- length(out$nocall); ptmp <- 100*ntmp/nts
	cat(ntmp," (",ptmp,"%) markers excluded because of low (<",callrate*100,"%) call rate\n",sep="")
	flush.console()
# check HWE
	if (!missing(hweidsubset)) {
		s <- summary(data[hweidsubset,])
	}
	Pexact <- s$Pexact
	rm(s);gc(verbose=FALSE)
	if (p.level < 0) {
		if (qoption == "storey") {
			if (require(qvalue)) {
				hweok <- !(qvalue(Pexact,fdr.level=fdrate)$significant)
			} else {
				warning("qvalue library not installed; using GenABEL's qvaluebh95() instead")
				hweok <- !(qvaluebh95(Pexact,fdrate)$significant)
			}
		} else {
			hweok <- !(qvaluebh95(Pexact,fdrate)$significant)
		}
	} else {
		hweok <- (Pexact >= p.level)
	}

# make output object
	out$nohwe <- data@snpnames[which(hweok == FALSE)]
#	cat("nohwe: ",out$nohwe,"\n")
	ntmp <- length(out$nohwe); ptmp <- 100*ntmp/nts
	if (p.level>=0)
		cat(ntmp," (",ptmp,"%) markers excluded because they are out of HWE (P <",p.level,")\n",sep="")
	else
		cat(ntmp," (",ptmp,"%) markers excluded because they are out of HWE (FDR <",fdrate,")\n",sep="")
	flush.console();
# all together
	out$snpok <- data@snpnames[(callok & freqok & redok & hweok)]
#	cat("ok: ",out$snpok,"\n")

# check call rate per person
	spp <- perid.summary(data[,out$snpok])
	idcallok <- (spp[,"CallPP"]>=perid.call)
	out$idnocall <- data@idnames[which(idcallok == FALSE)]
	ntmp <- length(out$idnocall); ptmp <- 100*ntmp/data@nids
	cat(ntmp," (",ptmp,"%) people excluded because of low (<",perid.call*100,"%) call rate\n",sep="")
	flush.console();
# check heterozygosity
	spp <- perid.summary(data[,out$snpok[out$snpok %in% autosomal(data)]])
	het <- spp[,"Het"]; 
	mh <- mean(het,na.rm=T); 
	sdh <- sd(het,na.rm=T)
	cat("Mean autosomal HET is ",mh," (s.e. ",sdh,")\n",sep="");
	flush.console()
	Zhet <- (het-mh)/sdh; 
	sigZ <- sign(Zhet)
	Zhet <- Zhet*Zhet
	Zhet <- (pchisq(Zhet,1,lower.tail=F))
	Zhet[is.na(Zhet)] <- 1
	if (qoption == "storey") {
		qobj <- qvalue(Zhet,fdr.level=(het.fdr*2.0))
		hetok <- (!(qobj$significant) | (sigZ<0 | is.na(sigZ)))
	} else {
		hetok <- (!(qvaluebh95(Zhet,(het.fdr*2.0))$significant) | (sigZ<0 | is.na(sigZ)))
	}
	if (any(!hetok)) {
		out$hetfail <- data@idnames[which(!hetok)]
		ntmp <- length(out$hetfail); ptmp <- 100*ntmp/data@nids
		cat(ntmp," (",ptmp,"%) people excluded because too high autosomal heterozygosity (FDR <",het.fdr*100,"%)\n",sep="")
		cat("Excluded people had HET >= ",min(het[!hetok]),"\n",sep="")
		flush.console();
	} else {
		cat("0 people excluded because too high autosomal heterozygosity (FDR <",het.fdr*100,"%)\n",sep="")
	}
	rm(spp);gc(verbose=FALSE)

# check IBS
	if (ibs.mrk<1) {
		out$ibsfail <- NULL;
		ibsok <- rep(TRUE,data@nids)
	} else {
#		fromset <- match(autosomal(data),out$snpok)
		fromset <- autosomal(data)[autosomal(data) %in% out$snpok]
		if (length(fromset) > ibs.mrk) {
			ibsset <- sort(sample(x=fromset,size=ibs.mrk,replace=FALSE))
			saibs <- ibs(data[,ibsset],weight="no")
		} else {
			saibs <- ibs(data[,fromset],weight="no")
			ibs.mrk <- length(fromset)
		}
		saibs[upper.tri(saibs,diag=T)] <- NA
		mibs <- mean(as.vector(saibs),na.rm=T)
		sdibs <- sd(as.vector(saibs),na.rm=T)
		if (any(saibs>=ibs.threshold,na.rm=T)) {
			ibsfailpairs <- crnames(dimnames(saibs),(saibs>=ibs.threshold))
			if (ibs.exclude == "both") {
				ibsfail <- unique(c(ibsfailpairs[,1],ibsfailpairs[,2]))
			} else {
				precll <- perid.summary(data[unique(ibsfailpairs[,1]),])
				preclr <- perid.summary(data[unique(ibsfailpairs[,2]),])
				cll <- precll[ibsfailpairs[,1],]$CallPP
				clr <- preclr[ibsfailpairs[,2],]$CallPP
				lgr <- (cll >= clr)
				ibsfail <- unique(c(ibsfailpairs[lgr,2],ibsfailpairs[!lgr,1]))
			}
			out$ibsfail <- ibsfail
			ibsok <- !(data@idnames %in% ibsfail)
			rm(ibsfailpairs,ibsfail);gc(verbose=FALSE)
		} else {
			out$ibsfail <- NULL;
			ibsok <- rep(TRUE,data@nids)
		}
		rm(saibs);gc(verbose=FALSE)
		cat("Mean IBS is ",mibs," (s.e. ",sdibs,"), as based on ",ibs.mrk," autosomal markers\n",sep="")
		flush.console();
		ntmp <- length(out$ibsfail); ptmp <- 100*ntmp/data@nids
		cat(ntmp," (",ptmp,"%) people excluded because of too high IBS (>=",ibs.threshold,")\n",sep="")
		flush.console()
	}

	out$idok <- data@idnames[(idcallok & hetok & ibsok)]

	ntmp <- length(out$snpok); ptmp <- 100*ntmp/nts
	cat("In total, ",ntmp," (",ptmp,"%) markers passed all criteria\n",sep="")
	flush.console();

	ntmp <- length(out$idok); ptmp <- 100*ntmp/data@nids
	cat("In total, ",ntmp," (",ptmp,"%) people passed all criteria\n",sep="")
	flush.console();

#Chi2 for the ones out of HWE, sorted
	out$Pex.nohwe <- Pexact[match(out$nohwe,data@snpnames)]
#	cat("Pex.nohwe: ",out$Pex.nohwe,"\n")
	idx <- sort(out$Pex.nohwe,decreasing=FALSE,ind=TRUE)$ix
	out$nohwe <- out$nohwe[idx]
#	cat("nohwe: ",out$nohwe,"\n")
	out$Pex.nohwe <- out$Pex.nohwe[idx]
#	cat("Pex.nohwe: ",out$Pex.nohwe,"\n")
# what was the call
	out$call$call <- match.call()
	out$call$name <- data@snpnames
	out$call$ids <- data@idnames
	out$call$map <- data@map
	out$call$chromosome <- data@chromosome
# output
	class(out) <- "check.marker"
	out
}

