## =======================================
## Extracts distinct states from sequences
## =======================================

STS_to_SPELL <- function(seqdata, id=NULL, pdata=NULL, birthdate=NULL, with.missing=TRUE) {

	if (!inherits(seqdata,"stslist")){
		stop("data is NOT a state sequence object, see seqdef function to create one")
	}
	nbseq <- nrow(seqdata)

	sl <- seqlength(seqdata)
	sltot <- sum(sl)

	void <- attr(seqdata, "void")
	statl <- attr(seqdata, "alphabet")
	nr <- attr(seqdata, "nr")
	if(is.null(id)){
		id <- rownames(seqdata)
		if(is.null(id)){
			id <- 1:nrow(seqdata)
		}
	} else if(length(id)==1 && nrow(seqdata)>1){
		if(is.null(pdata)){
			stop(" [!] no pdata provided and the 'id' argument seems to be a column name.")
		}
		id <- pdata[, id]
	}
	if(length(id)!=nrow(seqdata)){
		stop(" [!] 'id' should have one entry per sequence.")
	}
	if(!is.null(birthdate)){
		if(length(birthdate)==1 && nrow(seqdata)>1){
			if(is.null(pdata)){
				stop(" [!] no pdata provided and the 'birthdate' argument seems to be a column name.")
			}
			birthdate <- pdata[, birthdate]
		}
		
		if(length(birthdate)!=nrow(seqdata)){
			stop(" [!] 'birthdate' should have one entry per sequence.")
		}
		birthdate <- birthdate - 1
	}else{
		birthdate <- rep(0, nrow(seqdata))
	}
	begin <- numeric(sltot)
	end <-  numeric(sltot)
	ids <- vector(mode=mode(id), length=sltot)
	states <- character(sltot)
	if(with.missing) {
		statl <- c(statl, nr)
	}

	seqdatamat <- as.matrix(seqdata)
	
	if (!with.missing){
		seqdatamat[seqdatamat==nr] <- void
	}
	itot <- 1
	for (i in 1:nbseq) {
		
		idx <- 1
		sli <- sl[i]-1
		while (idx <= sl[i]) {
			ids[itot] <- id[i]
			iseq <- seqdatamat[i, idx]
			begin[itot] <- birthdate[i]+idx
			# if(itot ==1){
				# print(iseq)
				# print(str(states))
			# }
			while (idx <= sli && (seqdatamat[i, idx+1]==iseq)) {
				idx <- idx+1
			}

			if (iseq != void) {
				states[itot] <- as.character(iseq)
				end[itot] <- birthdate[i] + idx
				# if(itot ==1){
					# print(head(spell))
				# }
				itot <- itot+1
			}
			idx <- idx+1
		}

	}
	## drop=FALSE ensures that the result is a matrix even if trans has only one row
	keep <- 1:(itot-1)
	spell <- data.frame(id=ids[keep], begin=begin[keep], end=end[keep], states=factor(states[keep], levels=statl))

	return(spell)
}
