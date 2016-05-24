NestedDesign <- function(x, nlevel , indices = NULL, n = NULL){
	nx <- dim(as.matrix(x))[1]
	if(is.null(indices)){
		if(is.null(n)){
			warning("'indices' and 'n' cannot be both NULL")
		}else{
			indices <- list()
			indices[[1]] <- sort(sample(1:nx)[1:n[1]])
			if(nlevel > 2){
				for(i in 2:(nlevel-1)){
					indices[[i]] <- sort(sample(1:n[i-1])[1:n[i]])
				}
			}
		}
	}else{
		n <- c()
		for(i in 1:(nlevel-1)){
			n[i] <- length(indices[[i]])
		}	
	}
	MuFiDes <- list()
	MuFiDes$PX <- as.matrix(x)
	MuFiDes$ind <- indices
	MuFiDes$n <- n

	class(MuFiDes) <- "NestDesign"

	return(MuFiDes)

}

ExtractNestDesign <- function (NestDes, level) 
{
    ind <- NestDes$ind[[1]]
    if (level > 2) {
        for (l in 2:(level - 1)) {
            ind <- ind[NestDes$ind[[l]]]
        }
    }


    return(as.matrix(NestDes$PX[ind, ]))
}
