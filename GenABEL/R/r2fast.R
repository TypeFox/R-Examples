"r2fast" <- 
function(data,snpsubset,cross.snpsubset,idsubset) {
	if (is(data,"gwaa.data")) data <- data@gtdata;
	if (!is(data,"snp.data")) stop("The data argument must have snp.data-class");
	if (any(data@chromosome=="X") & dim(table(data@male))>1) {
		data <- data[,data@chromosome!="X"]
		warning("X-chromosome data dropped")
	}
###
	if (!missing(idsubset)) data <- data[idsubset,]
	if (missing(snpsubset) && !missing(cross.snpsubset)) stop("cross.snpsubset arg cannot be used (snpsubset missing)",immediate. = TRUE)
	r2.C.option <- 0
	if (!missing(snpsubset) && !missing(cross.snpsubset)) {
		snpset1 <- snpsubset
		snpset2 <- cross.snpsubset
		if (any(snpset1 %in% snpset2)) stop("snpsubset and cross.snpsubset should not overlap!")
		snpsorder <- c(snpset1,snpset2)
		if (length(snpsorder) != data@nsnps) data <- data[,snpsorder]
		if (any(snpsorder!=data@snpnames)) data <- data[,snpsorder]
		r2.C.option <- 1
	} else if (!missing(snpsubset) & missing(cross.snpsubset)) {
		snpset1 <- snpsubset
		snpset2 <- snpsubset
		snpsorder <- snpsubset
		data <- data[,snpsorder]
	} else if (missing(snpsubset) & missing(cross.snpsubset)) {
		snpset1 <- data@snpnames
		snpset2 <- data@snpnames
	} else {
		stop("can not be: impossible combination of snpsubset and cross.snpsubset")
	}
	gc()
	snpset1.num <- which(data@snpnames %in% snpset1)
	snpset2.num <- which(data@snpnames %in% snpset2)
###
	if (r2.C.option==0) {
		out <- .C("r2",as.raw(data@gtps),as.integer(data@nids),as.integer(data@nsnps),out=double(data@nsnps*data@nsnps))$out;
		dim(out) <- c(data@nsnps,data@nsnps)
		out <- t(out)
		diag(out) <- NA
		colnames(out) <- data@snpnames
		rownames(out) <- data@snpnames
	} else if (r2.C.option==1) {
		sout <- .C("r2new",as.raw(data@gtps),as.integer(data@nids),as.integer(data@nsnps),as.integer(length(snpset1.num)),as.integer(snpset1.num-1),as.integer(length(snpset2.num)),as.integer(snpset2.num-1),out=double(2*length(snpset1.num)*length(snpset2.num)))$out;
		out <- list()
		out$num <- sout[1:(length(snpset1.num)*length(snpset2.num))]
		out$r2 <- sout[(length(snpset1.num)*length(snpset2.num)+1):(length(snpset1.num)*length(snpset2.num)*2)]
		dim(out$num) <- c(length(snpset1.num),length(snpset2.num))
		dim(out$r2) <- c(length(snpset2.num),length(snpset1.num))
		out$r2 <- t(out$r2)
		out$num <- t(out$num)
		rownames(out$num) <- snpset2
		colnames(out$num) <- snpset1
		rownames(out$r2) <- snpset1
		colnames(out$r2) <- snpset2
	} else {
		stop("can not be: incorrect r2.C.option")
	}
	out
}
