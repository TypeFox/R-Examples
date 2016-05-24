## ================================
## Convert from SPELL to STS format
## ================================

SPELL_to_STS <- function(seqdata) {

	lid <- unique(seqdata[,1])
	nbseq <- length(lid)

	trans <- matrix("", nrow=nbseq, ncol=1)
		
	for (i in 1:nbseq) {
		## OBLIGE D'UTILISER L'INDEX POUR LES COLONNES SINON NE MARCHE PAS
		spell <- seqdata[seqdata[,1]==lid[i],]
		idxmax <- nrow(spell)
		
		if (idxmax>0) {
			for (j in 1:idxmax) {
				dur <- spell[j,4]-spell[j,3]
				if (dur<=0 | is.na(dur)) {
					warning("negative or invalid duration for id ",lid[i],", spell ",j)
					seq <- "NA"
				}
				else	seq <- paste(rep(spell[j,5],dur),collapse="-")

				if (j==1) trans[i] <- paste(trans[i],seq,sep="") 
				else trans[i] <- paste(trans[i],seq,sep="-")
			}
		}
	}

	return(trans)
}
