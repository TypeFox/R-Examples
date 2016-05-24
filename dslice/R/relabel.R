relabel <- function(x)
{
	na_id <- which(is.na(x))
	if(length(na_id) > 0){
		warning("There exists NA")
		ux <- sort(unique(x[!is.na(x)]))
		if(is.numeric(x)){
			ux <- c(sum(ux), ux)
			x[na_id] = ux[1]
		}
		if(is.character(x)){
			nastr <- NULL
			for(i in 1:length(ux)){
				nastr <- paste(nastr, ux[i], sep="")
			}
			ux <- c(nastr, ux)
			x[na_id] = ux[1]
		}
	}else{
		ux <- sort(unique(x))
	}
	y <- vector(length = length(x), mode = "integer")
	for(i in 1:length(x)){
		y[i] <- which(ux == x[i]) - 1L
	}
	return(y)
}
