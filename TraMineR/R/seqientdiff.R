## =======================
## Within Sequence Entropy
## =======================

seqientdiff <- function(seqdata, norm=TRUE) {

	if (!inherits(seqdata, "stslist"))
		stop("data is NOT a sequence object, see seqdef function to create one")

    entropydiff <- function(dur, norm) {
		len <- sum(dur)
		## If we don't check, the funcion returns NA (division by 0)
		ent <- entropy(dur)
		if(ent>0 && norm) {
			## The maximum entropy is when length of the DSS=length of the sequence 
			p <- 1/len
			entmax <- (-len)*(p*log(p))
			ent <- ent/entmax
		}
		return(ent)
    }
	iseqtab <- seqdur(seqdata)
	iseqtab[is.na(iseqtab)] <- 0
	ient <- apply(iseqtab,1,entropydiff, norm=norm)
    
	ient <- as.matrix(ient)
	colnames(ient) <- "Hdss"
	rownames(ient) <- paste("[",seq(1:length(ient)),"]",sep="")
	return(ient)
}	
