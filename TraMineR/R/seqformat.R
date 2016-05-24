## ============================================
## TRANSLATION BETWEEN SEQUENCE REPRESENTATIONS
## ============================================

seqformat <- function(data, var=NULL, id=NULL,
	from, to, compressed=FALSE,
	nrep=NULL, tevent, stsep=NULL, covar=NULL,
	SPS.in=list(xfix="()", sdsep=","),
	SPS.out=list(xfix="()", sdsep=","),
	begin=NULL, end=NULL, status=NULL,
	process=TRUE, pdata=NULL, pvar=NULL, limit=100, overwrite=TRUE,
	fillblanks=NULL, tmin=NULL, tmax=NULL, nr="*") {

	## Checking the format
	list.from <- c("STS","SPS","SPELL")
	list.to <- c("STS","SPS","DSS","SRS","TSE")

	if (inherits(data,"stslist")) {
		message(" [>] input is a sequence object, converting from STS format")
		from <- "STS"
	} else if (missing(from) || !from %in% list.from)	
		stop("input format must be one of: ", paste(list.from), call.=FALSE)

	if (missing(to) || !to %in% list.to)	
		stop("output format must be one of: ", paste(list.to), call.=FALSE)

	##
	if (!is.null(covar)) covariates <- subset(data,,covar)

	if (!is.null(id)) ident <- as.matrix(subset(data,,id))
	else ident <- NULL

	## ========
	## From STS
	## ========
	if (from=="STS") {
		## Extracting the sequences from the data set
		seqdata <- seqxtract(data, var)

		if (inherits(data,"stslist")) {
			message(" [>] converting special codes for missing states to NA's")
			seqdata[seqdata==attr(data,"nr")] <- NA
			seqdata[seqdata==attr(data,"void")] <- NA
		}
		else {
			if (is.null(stsep)) {
				stsep <- seqfcheck(seqdata)
				if (stsep %in% c("-",":")) seqdata <- seqdecomp(seqdata,sep=stsep)
			}
			else {
				seqdata <- seqdecomp(seqdata,sep=stsep)
			}
		}

		trans <- seqdata
	}

	## ========
	## From SPS
	## ========
	else if (from=="SPS") {
		## Extracting the sequences from the data set
		seqdata <- seqxtract(data, var)

		if (is.null(stsep))	{	
			stsep <- seqfcheck(seqdata)
			if (stsep %in% c("-",":")) seqdata <- seqdecomp(seqdata,sep=stsep)
		}
		else {
			seqdata <- seqdecomp(seqdata,sep=stsep)
		}

		trans <- SPS_to_STS(seqdata, spsformat=SPS.in, nr=nr)
	}

	## ==========
	## From SPELL
	## ==========
	else if (from=="SPELL") {
		if (!is.null(var)) data <- subset(data,,var)

		if (!is.null(id)) ident <- as.matrix(subset(data,,id))
		else ident <- as.matrix(subset(data,,1))

		if (!is.null(begin)) sp2 <- as.matrix(subset(data,,begin))
		else sp2 <- as.matrix(subset(data,,2))

		if (!is.null(end)) sp3 <- as.matrix(subset(data,,end))
		else sp3 <- as.matrix(subset(data,,3))

		if (!is.null(status)) sp4 <- as.matrix(subset(data,,status))
		else sp4 <- as.matrix(subset(data,,4))

		seqdata <- data.frame(ident, sp2, sp3, sp4)

		## Extracting the sequences from the data set
		## seqdata <- seqxtract(data, var, data.frame=TRUE)
		
		trans <- BIOSPELL_to_STS(seqdata=seqdata,
			process=process, pdata=pdata, pvar=pvar,
			limit=limit, overwrite=overwrite, fillblanks=fillblanks,
			tmin=tmin, tmax=tmax)	
		## ident <- unique(ident)
	}

	## ===============
	## INTERNAL FORMAT
	## ===============
	rm(seqdata)
	nbin <- nrow(trans)
	if (from != "STS") message(" [>] ", from," data converted into ",nbin," STS sequences")

	## =============
	## OUTPUT FORMAT
	## =============
	if (to=="SPS") {
		out <- STS_to_SPS(seqdata=trans, spsformat=SPS.out)
		nbout <- seqdim(out)[1]

		if (compressed) out <- seqconc(out)
		}

	## To Distinct-State-Sequence format
	else if (to=="DSS") {
		out <- STS_to_DSS(trans)
		nbout <- seqdim(out)[1]

		if (compressed) out <- seqconc(out)
	}

	## STS
	else if (to=="STS") {
		out <- trans
		nbout <- nrow(out)

		if (compressed) out <- seqconc(out)
		}

	else if (to=="SRS") {
		out <- STS_to_SRS(trans,nrep)
		if (!is.null(covar)) out <- merge(out,data.frame(id=seq(1:nbin),covariates))
		nbout <- nrow(out)
	}

	else if (to=="TSE") {
		if (missing(tevent))
			stop(" [!] you must provide a transition-definition matrix (see seqetm)", call.=FALSE)
		out <- STS_to_TSE(trans, ident, tevent)
		nbout <- nrow(out)
		}
	else
		stop("Output format is unsupported")

	if (to!="STS") message(" [>] STS sequences converted to ",nbout," ", to," seq./rows")
	return(out)

}
