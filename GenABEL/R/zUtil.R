# utility functions invisible to user
lossFunctionLambdaKS <- function(lambda,chi2values, ... ) {
	ksT <- ks.test(chi2values/lambda, ... )$stat
	return( ksT )
}
estLambdaKS <- function(chi2values,limits=c(0.5,100),df=1) {
	iniLambda <- 1
	optRes <- optimize(lossFunctionLambdaKS, interval=limits, chi2values=chi2values, "pchisq", 1,df=df)
	return(optRes$minimum)
}
chrom.char2num <- function(chrom) {
	chrom <- as.character(chrom)
	chdesU <- unique(chrom)
	chnumU <- suppressWarnings(as.numeric(chdesU))
	chCHU <- sort(chdesU[is.na(chnumU)])
	if (any(!is.na(chnumU))) {
		maxch <- max(chnumU,na.rm=T)
	} else {
		maxch <- 0
	}
	j <- 1
	for (i in chCHU) {
		chrom[which(chrom==i)] <- as.character(maxch+j)
		j <- j + 1
	}
	chnum <- as.numeric(chrom)
	chnum
}
