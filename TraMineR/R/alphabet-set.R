## ====================================
## Function to set the alphabet
## of a sequence object
## ==================================== 

"alphabet<-" <- function(seqdata, value) {

	if (!inherits(seqdata,"stslist")) 
		stop("data is not a sequence object, see seqdef function to create one")

	nbstate <- length(alphabet(seqdata))
	
	if (length(value)!=nbstate) 
		stop("number of states must be ",nbstate)
	else {
		for (i in 1:dim(seqdata)[2])
			levels(seqdata[,i])[1:nbstate] <- value

		if (all.equal(attr(seqdata,"alphabet"),attr(seqdata,"labels")))
			attr(seqdata,"labels") <- value

		attr(seqdata,"alphabet") <- value
	}

	seqdata

} 
