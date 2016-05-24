"redundant" <-
function(data,pairs = "bychrom", minconcordance = 2.0) {
#data is snp.data
	if (!is(data,"snp.data")) stop("Wrong class of data (should be snp.data)");
	lst <- new.env()
	if (pairs == "all") {
		lst[["all"]] <- c(1:data@nsnps)
		clev <- "all"
	} else if (pairs == "bychrom") {
		clev <- levels(as.factor(data@chromosome))
		for (i in 1:length(clev)) 
			lst[[clev[i]]] <- which(data@chromosome == clev[i])
	} else {
		stop("pairs argument should be \"all\" or \"bychrom\"")
	}
	lst <- as.list(lst)
	out <- new.env()
	for (i in 1:length(lst)) {
		n1data <- data[,lst[[clev[i]]]]
		idx <- .C("redundant",as.raw(n1data@gtps),as.integer(n1data@nids),as.integer(n1data@nsnps), as.double(minconcordance), outlist = integer(n1data@nsnps), PACKAGE="GenABEL")$outlist
		red <- which(idx>0)
		red <- n1data@snpnames[red]
		out[["all"]] <- c(out[["all"]],red)
		for (i in 1:length(idx)) {
			if (idx[i]) {
				nam <- n1data@snpnames[idx[i]]
				out[[nam]] <- c(out[[nam]],n1data@snpnames[i])
			}
		}
		out <- as.list(out)
	}
	out
}

