## SEQUENCE TRUNCATION

TraMineR.trunc <- function(seq, mstate, sl, 
	left="DEL", right="DEL", gaps="DEL", 
	neutral="#", void="%") {

	sidx <- 1:sl

	## Index des missing et index des etats valides
	na.pos <- sidx[mstate]
	notna.pos <- sidx[!mstate]
	
	## Position du premier etat valide
	c1 <- notna.pos[1]

	if (c1>1)
		lc <- c1-1
	else lc=0

	rc <- max(notna.pos)+1
	mm <- na.pos[na.pos > lc+1 & na.pos < rc-1]

	seq.trunc <- seq
	
	if (!is.na(left) & lc>0) {
		if (left=="DEL") seq.trunc[1:lc] <- void
		else if (left=="NEUTRAL") seq.trunc[1:lc] <- neutral
		else seq.trunc[1:lc] <- left
	}

	if (!is.na(right) & rc<=sl) {
		if (right=="DEL") seq.trunc[rc:sl] <- void
		else if (right=="NEUTRAL") seq.trunc[rc:sl] <- neutral
		else seq.trunc[rc:sl] <- right
	}

	if (!is.na(gaps) & length(mm>0)) {
		if (gaps=="DEL") seq.trunc[mm] <- void
		else if (gaps=="NEUTRAL") seq.trunc[mm] <- neutral
		else seq.trunc[mm] <- gaps
	}
	
	ndel <- sum(seq.trunc==void, na.rm=TRUE)
	
	if (ndel>0) {
		seq.trunc <- seq.trunc[seq.trunc!=void]
		seq.trunc <- c(seq.trunc,rep(void,ndel))
	}

	return(seq.trunc)
} 
