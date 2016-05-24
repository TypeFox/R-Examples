binarizeMatrix <- function(mat, method=c("BASCA","BASCB","kMeans"), adjustment="none", ...){
	binFunc <- switch(match.arg(method, c("BASCA","BASCB","kMeans")),
		"BASCA" = function(x, ...){
			bin <- binarize.BASC(x, method="A", ...)
			return(c(bin@binarizedMeasurements, 
				list(bin@threshold, bin@p.value)))
		},
		"BASCB" = function(x, ...){
			bin <- binarize.BASC(x, method="B", ...)
			return(c(bin@binarizedMeasurements, 
				list(bin@threshold, bin@p.value)))
		},
		"kMeans" = function(x, ...){
			bin <- binarize.kMeans(x, ...)
			return(c(bin@binarizedMeasurements, 
				list(bin@threshold, bin@p.value)))
		})

	bin <- do.call("rbind.data.frame", apply(mat, 1, binFunc, ...))

	if(!is.null(colnames(mat))){
		colnames(bin) <- c(colnames(mat), "threshold", "p.value")
	}else{
		colnames(bin) <- c(paste("V",seq_len(ncol(mat)),sep=""),
						"threshold", "p.value")
	}

	bin[,"p.value"] <- p.adjust(bin[,"p.value"], method=adjustment)
	return(bin)
}

trinarizeMatrix <- function(mat, method=c("TASCA","TASCB","kMeans"), adjustment="none", ...){
	triFunc <- switch(match.arg(method, c("TASCA","TASCB","kMeans")),
		"TASCA" = function(x, ...){
			tri <- TASC(x, method="A", ...)
			return(c(tri@trinarizedMeasurements, 
				list(tri@threshold1, tri@threshold2, tri@p.value)))
		},
		"TASCB" = function(x, ...){
			tri <- TASC(x, method="B", ...)
			return(c(tri@trinarizedMeasurements, 
				list(tri@threshold1, tri@threshold2, tri@p.value)))
		},
		"kMeans" = function(x, ...){
			tri <- trinarize.kMeans(x, ...)
			return(c(tri@trinarizedMeasurements, 
				list(tri@threshold1, tri@threshold2, tri@p.value)))
		})

	tri <- do.call("rbind.data.frame", apply(mat, 1, triFunc, ...))

	if(!is.null(colnames(mat))){
		colnames(tri) <- c(colnames(mat), "threshold1", "threshold2", "p.value")
	}else{
		colnames(tri) <- c(paste("V",seq_len(ncol(mat)),sep=""),
						"threshold1", "threshold2", "p.value")
	}

	tri[,"p.value"] <- p.adjust(tri[,"p.value"], method=adjustment)
	return(tri)
}
