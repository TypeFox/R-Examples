"snp.data" <-
function (nids,rawdata,
		idnames=as.character(c(1:nids)),
		snpnames=as.character(c(1:(length(rawdata)/ceiling(nids/4)))),
		chromosome=as.factor(c(1:(length(rawdata)/ceiling(nids/4)))),
		map=as.double(seq(1,(length(rawdata)/ceiling(nids/4)))),
		coding=as.raw(rep(1,length(rawdata)/ceiling(nids/4))),
		strand=as.raw(rep(0,length(rawdata)/ceiling(nids/4))),
		male=rep(0,nids)
		) {
	nbytes <- ceiling(nids/4)
	nsnps <- length(rawdata)/nbytes
	if (nsnps != round(nsnps)) stop("wrong numner of ids")
	if (class(rawdata) != " snp.mx") stop("rawdata must have class snp.mx")
#	t <- as.raw(rawdata)
#	dim(t) <- c(nbytes,nsnps)
#	t <- new("snp.mx",t)
	a <- new("snp.data",nbytes=nbytes,nids=nids,nsnps=nsnps,gtps=rawdata,
			idnames=idnames,snpnames=snpnames,
			chromosome=chromosome,map=map,coding=coding,
			strand=strand,male=male)
	a
}

