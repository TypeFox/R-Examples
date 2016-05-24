## ============================================
## Returns the alphabet for a sequence data set
## ============================================

seqstatl <- function(data, var=NULL, format='STS') {

	listform <- c("STS","SPS","SPELL","DSS")
	if (!format %in% listform) 
		stop("Input format must be one of: ", listform)

	## Extracting the sequences from the data set
	seqdata <- seqxtract(data, var)

	## 
	if (format=='SPS') seqdata <- seqformat(seqdata, from='SPS', to='STS')

	## Convert into the extended format to list states/events
	if (seqfcheck(seqdata) %in% c(":","-"))
		seqdata <- seqdecomp(seqdata, sep=seqfcheck(seqdata))

	statl <- sort(unique(as.vector(seqdata)))

	## IF states are numeric values, sort them as integer
	## (if sorted as characters, 10 will be after 1,
	## not after 9)
	statnum <- suppressWarnings(as.integer(statl))

	if (all(is.na(statnum)==FALSE)) statl <- sort(statnum)

	return(statl)
	}

