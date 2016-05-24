## ================================================
## Retrieve the state labels from a sequence object
## ================================================

stlab <- function(seqdata) {
	statelab <- attr(seqdata,"labels")
	return(statelab)
}

