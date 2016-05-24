"ccfast.new" <-
function(y,data,snpsubset,idsubset,quiet=FALSE) {
	if (!is(data,"gwaa.data")) stop("wrong type of data argument, must be gwaa.data")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (any(data@gtdata@chromosome=="X") & dim(table(data@gtdata@male))>1) {
		data <- data[,data@gtdata@chromosome!="X"]
		if (!quiet) cat("X-chromosome data dropped\n")
	}

#	attach(data@phdata,warn.conflicts=FALSE,pos=2)
#	cc <- get(y,pos=2)
#	detach(data@phdata)
	cc <- phdata(data)[[y]]

        if (length(levels(as.factor(cc)))<2) stop("cc status is monomorphic!") 
        if (length(levels(as.factor(cc)))>2) stop("cc status has more than 2 levels!") 
        if (levels(as.factor(cc))[1] != 0 || levels(as.factor(cc))[2] != 1) stop ("cc is case-control status, with 0 as control and 1 as cases. No 0 and/or 1 found in the data")

	if (any(is.na(cc))) {
	  if (!quiet) warning(paste(sum(is.na(cc)),"people (out of",length(cc),") excluded as not having cc status\n"),immediate. = TRUE)
	  vec <- !is.na(cc)
	  data <- data[vec,]
	  cc1 <- cc[!is.na(cc)]
	} else {
	  cc1 <- cc
	}
	rm(cc)

	a <- fcc.new(data@gtdata,cc1)
        lena <- data@gtdata@nsnps

	out <- list()
	out$Padd <- a[1:lena] #(1. - pchisq(a[1:lena],1))
#	out$medadd <- median(a[1:lena])
	out$Pdom <- a[(lena+1):(2*lena)] #(1. - pchisq(a[(lena+1):(2*lena)],1))
#	out$meddom <-  median(a[(lena+1):(2*lena)])
	out$Prec <- a[(2*lena+1):(3*lena)] #(1. - pchisq(a[(2*lena+1):(3*lena)],1))
#	out$medrec <-  median(a[(2*lena+1):(3*lena)])
	out$effadd <- a[(3*lena+1):(lena*4)]
	out$effdom <- a[(4*lena+1):(lena*5)]
	out$effrec <- a[(5*lena+1):(lena*6)]
	rm(a);gc(verbose=FALSE)
	out$name <- data@gtdata@snpnames
	out$formula <- match.call()
	out$family <- "chi-square 1 and 2 d.f."
	out$map <- data@gtdata@map
	out$chromosome <- data@gtdata@chromosome
	out$ids <- data@gtdata@idnames
#	out$lambda <- estlambda(out$P1df,plot=FALSE,prop=0.9)
	class(out) <- "scan.gwaa"
	out
}

