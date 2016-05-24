## ====================================
## Function to set the state labels
## of a sequence object
## ==================================== 

"stlab<-" <- function(seqdata, value) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one")

	nbstate <- length(alphabet(seqdata))
	
	if (length(value)!=nbstate) 
		stop("number of labels must be",nbstate)
	else
		attr(seqdata,"labels") <- value

	seqdata
} 
