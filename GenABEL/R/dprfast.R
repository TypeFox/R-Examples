"dprfast" <- 
function(data,snpsubset,idsubset) {
	if (is(data,"gwaa.data")) data <- data@gtdata;
	if (!is(data,"snp.data")) stop("The data argument must have snp.data-class");
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (any(data@chromosome=="X") & dim(table(data@male))>1) {
		data <- data[,data@chromosome!="X"]
		warning("X-chromosome data dropped")
	}
	out <- .C("dprime",as.raw(data@gtps),as.integer(data@nids),as.integer(data@nsnps),out=double(data@nsnps*data@nsnps))$out;
	dim(out) <- c(data@nsnps,data@nsnps)
	out <- t(out)
	diag(out) <- NA
	colnames(out) <- data@snpnames
	rownames(out) <- data@snpnames
	out
}
