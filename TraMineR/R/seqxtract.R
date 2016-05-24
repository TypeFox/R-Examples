## =================================
## EXTRACT SEQUENCES FROM A DATA SET
## =================================

seqxtract <- function(data, var, data.frame=FALSE) {
	
	## Extracting the sequences from the data set
	if (missing(var) || is.null(var) || is.na(var)) 
		seqdata <- data
	else 
		seqdata <- subset(data,,var)

	if (data.frame==FALSE) 
		seqdata <- as.matrix(seqdata)

	return(seqdata)

	}
