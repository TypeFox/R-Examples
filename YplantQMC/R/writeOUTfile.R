writeOUTfile <- function(yplantsim){

	if(!inherits(yplantsim, "yplantsim"))
		stop("Need object output by YplantDay.")
	
	pfile <- yplantsim$plant$pfile
	proot <- gsub("\\.p$","",pfile, ignore.case=TRUE)
	filen <- paste0(proot,"-YplantQMC.OUT")
	
	write.table(yplantsim$outdata, filen, sep="\t", quote=FALSE, row.names=FALSE)
	message("Results written to file: ", filen)

return(invisible(filen))
}