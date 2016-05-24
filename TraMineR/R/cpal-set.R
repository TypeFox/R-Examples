## ====================================
## Function to set the color palette
## of a sequence object
## ==================================== 

"cpal<-" <- function(seqdata, value) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one")

	nbstate <- length(alphabet(seqdata))
	
	if (length(value)!=nbstate) 
		stop("number of colors must be",nbstate)
	else
		attr(seqdata,"cpal") <- value

	seqdata
} 
