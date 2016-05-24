"Xfix" <-
function(data) {
	if (is(data,"gwaa.data")) data@gtdata <- Xfix.internal(data@gtdata)
	else if (is(data,"snp.data")) data <- Xfix.internal(data)
	else stop("data argument must be of gwaa.data-class or snp.data-class")
	data
}

"Xfix.internal" <-
function(data) {
	if (!is(data,"snp.data")) stop("data argument must be of snp.data-class")
	if (!any(data@chromosome == "X") & !any(data@chromosome == "Y") & !any(data@chromosome == "mt")) stop("No X-, Y- or mtDNA-linked markers")
#	mlst <- data@male==1
#	xmrk <- data@chromosome=="X"
#	if (sum(mlst)<1 || sum(xmrk)<1) {
#		cat("no SNPs typed in male\n")
#		return(data)
#	}
#	print("entering imphetckeck")
	err <- imphetcheck(data)
#	print("exiting imphetckeck")
#	print(dim(err$Xerrtab))
#	print(err)
	if (!err$xerr) {
		cat("no X/Y/mtDNA-errors to fix\n")
		return(data)
	} else {
#		print(dim(err$Xerrtab))
		fsnps <- unique(err$Xerrtab[,2])
#		print(length(fsnps))
		errlist <- list()
		errlist[["X"]] <- 0
		errlist[["Y"]] <- 0
		errlist[["mt"]] <- 0
		for (i in fsnps) {
#			print(i)
			fids <- unique(err$Xerrtab[err$Xerrtab[,2]==i,1])
			gtv <- as.numeric(data[,i])+1
			gtv[is.na(gtv)] <- 0
			gtv[(data@idnames %in% fids)] <- 0
#			print(dim(data@gtps[,which(data@snpnames==i)]))
#			print(dim(gtv))
#			print(length(data@gtps[,which(data@snpnames==i)]))
#			print(length(gtv))
			errlist[[as.character(data@chromosome[i])]] <- errlist[[as.character(data@chromosome[i])]] + length(fids)
#			print(i)
#			print(as.character(data@chromosome[i]))
#			print(errlist)
			data@gtps[,which(data@snpnames==i)] <- put.snps(gtv)
		}
#		for (i in 1:dim(err$Xerrtab)[1]) {
#			cid <- which(data@idnames==err$Xerrtab[i,1])
#			csn <- which(data@snpnames==err$Xerrtab[i,2])
#			gtv <- as.numeric(data[,csn])+1
#			gtv[is.na(gtv)] <- 0
#			gtv[cid] <- 0
#			data@gtps[,csn] <- put.snps(gtv)
#		} 
		cat("...",(errlist$X+errlist$Y+errlist$mt),"X/Y/mtDNA (",errlist$X,errlist$Y,errlist$mt,") impossible heterozygotes and female Ys set as missing\n")
	}
	data
}
