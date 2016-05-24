print.select <- function(x, ...){

	## Print a summary of atom selection object features
  	if(!inherits(x, "select")) {
    	stop("Input should be a 'select' object, as obtained from 'atom.select()'")
  	}

 	cat("\n Call:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
 		"\n", sep = "")	

	cat( paste0("\n   Atom Indices#: ", length(x$atom), "  ($atom)",
		"\n   XYZ  Indices#: ",length(x$xyz), "  ($xyz)\n\n" ) )
	## Add call
	if(	length(x$atom) != (length(x$xyz)/3) ) 
		warning("Atom and XYZ Indices Miss-Match")

	i <- paste(attributes(x)$names, collapse = ", ")
	cat(strwrap(paste(" + attr:", i), width = 45, exdent = 8), sep = "\n")

}

