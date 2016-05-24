#'@method print plant3d
#'@S3method print plant3d
print.plant3d <- function(x, hint=TRUE, ...){

	cat("Yplant - plant object (class 'plant3d')\n\n")
	cat(paste(c(rep("-",30),"\n"),collapse=""))
	cat("Constructed from", x$inputformat, "file format\n")
	if(x$inputformat == "P")cat("Input files used :", x$pfile, x$lfile, "\n\n")
	if(x$inputformat == "Q"){
		if(is.character(x$qdata))cat("Input file used :", x$qdata, "\n")
		if(!is.null(x$qfile) && is.data.frame(x$qfile))
      cat("Input files used :\n", x$qfile,"\n", x$lfile, "\n")
	}
  cat("\n")
  
	if(x$inputformat == "P")cat("Plant has", x$nleaves, "leaves on", nrow(x$pdata), "nodes.\n")
	if(x$inputformat == "Q")cat("Plant has", x$nleaves, "leaves.\n")
	
	if(hint)cat("To generate a detailed summary of your plant, use the summary() function.\n")
}