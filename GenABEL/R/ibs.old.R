"ibs.old" <- 
function (data,snpsubset,idsubset,weight="no") {
	if (is(data,"gwaa.data")) data <- data@gtdata
	if (!is(data,"snp.data")) stop("The data argument must be of snp.data-class or gwaa.data-class")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	homodiag <- hom(data)
	if (!missing(idsubset)) homodiag <- homodiag[idsubset,]
	if (weight=="no") {
		homodiag <- homodiag[,"Hom"]
	} else {
		homodiag <- 0.5+homodiag[,"F"]
	}
	if (!missing(idsubset)) data <- data[idsubset,]
	wargs <- c("no","freq")
	if (!(match(weight,wargs,nomatch=0)>0)) {
		out <- paste("weight argument should be one of",wargs,"\n")
		stop(out)
	}
#	if (npairs > 2500000) stop("Too many pairs to evaluate... Stopped")
	if (weight == "no") option = 0
	if (weight == "freq") option = 1
	out <- .C("ibs", as.raw(data@gtps), as.integer(data@nids), as.integer(data@nsnps), as.integer(option), sout = double(data@nids*data@nids), PACKAGE="GenABEL")$sout
	dim(out) <- c(data@nids,data@nids)
	out <- t(out)
	diag(out) <- homodiag
	colnames(out) <- data@idnames
	rownames(out) <- data@idnames
	out
}
