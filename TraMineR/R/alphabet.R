## ============================================
## Retrieve the alphabet from a sequence object
## ============================================

alphabet <- function(seqdata) {
	statl <- attr(seqdata,"alphabet")
	return(statl)
}
	
