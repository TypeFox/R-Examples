"save.snp.data" <-
function(data, genofile = "geno.raw",human=FALSE) {
	if (is(data,"gwaa.data")) data <- data@gtdata
	if (!is(data,"snp.data")) stop("data argument must be of snp.data-class or gwaa.data-class!")
	if (human==TRUE && missing(genofile)) genofile="geno.dat"
	ofile <- file(genofile,"w")
	if (human) 
		cat("#GenABEL human-readable data version 0.1\n",file=ofile)
	else 
		cat("#GenABEL raw data version 0.1\n",file=ofile)
	cat("writing id names...\n")
	cat(data@idnames,"\n",file=ofile)
	cat("writing snp names...\n")
	cat(data@snpnames,"\n",file=ofile)
	cat("writing chromosome data...\n")
	cat(as.character(data@chromosome),"\n",file=ofile)
	cat("writing map data...\n")
	cat(data@map,"\n",file=ofile)
	cat("writing coding data...\n")
	if (human) 
		cat(as.character(data@coding),"\n",file=ofile)
	else 
		write(data@coding,file=ofile,ncolumns=data@nsnps)
	cat("writing strand data...\n")
	if (human) 
		cat(as.character(data@strand),"\n",file=ofile)
	else 
		write(data@strand,file=ofile,ncolumns=data@nsnps)
	cat("writing genotypic data...\n")
	if (human) {
		idta <- as.double(data)+1
		idta <- replace(idta,is.na(idta),0)
		write(idta,file=ofile,ncolumns=data@nids)
	} else 
		write(data@gtps,file=ofile,ncolumns=data@nbytes)
#	for (i in 1:data@nsnps) {
#		cat(data@gtps[,i],"\n",file=ofile)
#	}
	close(ofile)
}

