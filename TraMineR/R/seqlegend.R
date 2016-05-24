## =========================================
## Plotting the legend for a sequence object
## =========================================

seqlegend <- function(seqdata, with.missing="auto", cpal=NULL, missing.color=NULL, ltext=NULL, 
	position="topleft", fontsize=1,...) {
	
	if (!inherits(seqdata,"stslist"))
		stop("data is not a sequence object, use seqdef function to create one")

	if (is.null(cpal)) 
		cpal <- attr(seqdata,"cpal")

	if (is.null(ltext)) 
		ltext <- attr(seqdata,"labels")

	if (is.null(missing.color)) 
		missing.color <- attr(seqdata,"missing.color") 

	## Adding an entry for missing in the legend
	nr <- attr(seqdata,"nr")

	if ((with.missing=="auto" && any(seqdata==nr)) || with.missing==TRUE) {
		cpal <- c(cpal,missing.color)
		ltext <- c(ltext,"missing")
	## statl <- c(statl,nr)
	## nbstat <- nbstat+1
	}

	plot(0, type= "n", axes=FALSE, xlab="", ylab="")
	legend(position, fill=cpal, legend=ltext, cex=fontsize,...)
}
