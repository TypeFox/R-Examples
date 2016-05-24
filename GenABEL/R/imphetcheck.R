"imphetcheck" <-
function(data) {
	if (!is(data,"snp.data")) stop("data argument should be of snp.data-class")
	if (!any(data@chromosome == "X") & !any(data@chromosome == "Y") & !any(data@chromosome == "mt")) stop("No X-, Y- or mtDNA-linked markers")
	male <- (data@male==1)
	out <- list()
	out$xerr <- 0
# male X checks
	if (sum(male) & any(data@chromosome == "X")) {
		xdat <- as.numeric(data[male,(data@chromosome == "X")])
		if (any(xdat==1,na.rm=T)) {
			out$xerr <- 1
			out$Xerrtab <- crnames(dimnames(xdat),which(xdat==1))
		}
	}
# Y checks
	if (any(data@chromosome == "Y")) {
		xdat <- as.numeric(data[,(data@chromosome == "Y")])
		if (any(xdat==1,na.rm=T)) {
			if (!out$xerr) {
				out$xerr <- 1
				out$Xerrtab <- crnames(dimnames(xdat),which(xdat==1))
			} else {
				out$Xerrtab <- rbind(out$Xerrtab,crnames(dimnames(xdat),which(xdat==1)))
			}
		}
# female Y checks
		if (any(data@male==0)) {
		xdat <- as.numeric(data[(data@male==0),(data@chromosome == "Y")])
		if (any(!is.na(xdat))) {
			if (!out$xerr) {
				out$xerr <- 1
				out$Xerrtab <- crnames(dimnames(xdat),which(!is.na(xdat)))
			} else {
				out$Xerrtab <- rbind(out$Xerrtab,crnames(dimnames(xdat),which(!is.na(xdat))))
			}
		}
		}
	}
# mtDNA checks
	if (any(data@chromosome == "mt")) {
		xdat <- as.numeric(data[,(data@chromosome == "mt")])
		if (any(xdat==1,na.rm=T)) {
			if (!out$xerr) {
				out$xerr <- 1
				out$Xerrtab <- crnames(dimnames(xdat),which(xdat==1))
			} else {
				out$Xerrtab <- rbind(out$Xerrtab,crnames(dimnames(xdat),which(xdat==1)))
			}
		}
	}
# output
	if (out$xerr) colnames(out$Xerrtab) <- c("ID","SNP")
	out
}
