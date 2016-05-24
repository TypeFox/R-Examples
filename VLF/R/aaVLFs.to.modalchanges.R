aaVLFs.to.modalchanges <-
function(modal, AminoAcidList, aalength){
	polarcharged <- c("D", "E", "K", "R", "H")
	polaruncharged <- c("S", "T", "Q", "N", "Y")
	Nonpolar <- c("A", "V", "L", "I", "M", "F", "W")
	Unique <- c("G", "C", "P")

	positionAllAminos <- list()
	p <- 1
	for(i in 1:nrow(AminoAcidList)){
		for(n in 1:aalength){
			if(is.na(AminoAcidList[i, n+2]) == FALSE){
				positionAllAminos[[p]] <- c(AminoAcidList[i,1], AminoAcidList[i, 2], n, AminoAcidList[i,n+2])
				p = p + 1
			}
		}
	}

	sameAll <- 0
	changedAll <- 0
	z <- c()
	r <- c()
	for(i in 1:length(positionAllAminos)){
		z[i] <- (modal[as.numeric(positionAllAminos[[i]][3])])
		r[i] <- positionAllAminos[[i]][4]
	
		if(any(z[i] == polaruncharged) && any(r[i] == polaruncharged)){
			sameAll = sameAll + 1
			}else{
			if(any(z[i] == polarcharged) && any(r[i] == polarcharged)){
				sameAll = sameAll + 1
				}else{
				if(any(z[i] == Nonpolar) && any(r[i] == Nonpolar)){
					sameAll = sameAll + 1
					}else{
					if(any(z[i] == Unique) && any(r[i] == Unique)){
						sameAll = sameAll + 1
					}else{changedAll = changedAll + 1}
				}
			}	
		}
	}
	return(rbind(sameAll,changedAll))
}
