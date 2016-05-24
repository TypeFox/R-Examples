concordant.to.modalchanges <-
function(matched, modal){
	polarcharged <- c("D", "E", "K", "R", "H")
	polaruncharged <- c("S", "T", "Q", "N", "Y")
	Nonpolar <- c("A", "V", "L", "I", "M", "F", "W")
	Unique <- c("G", "C", "P")

	x<- c()
	y <- c()
	same = 0
	changed = 0
	continue = TRUE
	for(i in 1:length(matched)){
		if(i > 1 && matched[[i]][1] == matched[[i-1]][1]){
			if(matched[[i]][4] == matched[[i-1]][4]){
				continue = FALSE
			}else{continue = TRUE}
		}else{continue = TRUE}

		if(continue == TRUE){
			x[i] <- (modal[as.numeric(matched[[i]][4])])
			y[i] <- matched[[i]][3]
		}
		
		if(is.na(x[i]) == FALSE && any(x[i] == polaruncharged) && any(y[i] == polaruncharged)){
			same = same + 1
			}else{
			if(is.na(x[i]) == FALSE && any(x[i] == polarcharged) && any(y[i] == polarcharged)){
				same = same + 1
				}else{
				if(is.na(x[i]) == FALSE && any(x[i] == Nonpolar) && any(y[i] == Nonpolar)){
					same = same + 1
					}else{
					if(is.na(x[i]) == FALSE && any(x[i] == Unique) && any(y[i] == Unique)){
						same = same + 1
					}else{
					if(is.na(x[i]) == FALSE){
						changed = changed + 1
						}
					}
				}
			}	
		}
	}
	return(rbind(same, changed))
}
