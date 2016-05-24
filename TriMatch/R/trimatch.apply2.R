#' Recursive function to find possible matched triplets using the apply functions.
#' 
#' Internal method. This version does not use the exact matching. Instead, this
#' function should be called separately for each grouping.
#' 
#' @inheritParams trimatch
#' @param sd1 standard deviation for propensity scores from model 1.
#' @param sd2 standard deviation for propensity scores from model 2.
#' @param sd3 standard deviation for propensity scores from model 3.
trimatch.apply2 <- function(tpsa, caliper, nmatch, match.order, sd1, sd2, sd3) {	
	ps1 <- getPS(tpsa, match.order[1], match.order[2])
	ps2 <- getPS(tpsa, match.order[2], match.order[3])
	ps3 <- getPS(tpsa, match.order[3], match.order[1])
	
	#Calculate the distances
	ps1.2 <- ps1[tpsa$treat == match.order[2]]
	ps1.1 <- ps1[tpsa$treat == match.order[1]]
	d1 <- abs(sapply(ps1.2, FUN=function(x) { x - ps1.1 })) / sd1		
	dimnames(d1) <- list(
		tpsa[tpsa$treat == match.order[1],'id'],
		tpsa[tpsa$treat == match.order[2],'id']	)
	
	ps2.3 <- ps2[tpsa$treat == match.order[3]]
	ps2.2 <- ps2[tpsa$treat == match.order[2]]
	d2 <- abs(sapply(ps2.3,FUN=function(x) { x - ps2.2 })) / sd2
	dimnames(d2) <- list(
		tpsa[tpsa$treat == match.order[2],'id'],
		tpsa[tpsa$treat == match.order[3],'id']	)
	
	ps3.1 <- ps3[tpsa$treat == match.order[1]]
	ps3.3 <- ps3[tpsa$treat == match.order[3]]
	d3 <- abs(sapply(ps3.1, FUN=function(x) { x - ps3.3 })) / sd3
	dimnames(d3) <- list(
		tpsa[tpsa$treat == match.order[3],'id'],
		tpsa[tpsa$treat == match.order[1],'id']	)

	#Check the caliper
	d1[d1 > caliper[1]] <- NA
	d2[d2 > caliper[2]] <- NA
	d3[d3 > caliper[3]] <- NA
	
	match1 <- function(i1) {
		row1 <- d1[i1,]
		row1 <- row1[!is.na(row1)]
		row1 <- row1[order(row1)]
		nxt <- names(row1)[seq_len(min(nmatch[1], length(row1)))]
		if(length(nxt) > 0) {
			return(lapply(nxt, match2, i1=i1, row1=row1))
		} else {
			return(c())
		}
	}
	
	match2 <- function(i2, i1, row1) {
		row2 <- d2[i2,]
		row2 <- row2[!is.na(row2)]
		row2 <- row2[order(row2)]
		nxt <- names(row2)[seq_len(min(nmatch[2], length(row2)))]
		if(length(nxt) > 0) {
			return(lapply(nxt, match3, i1=i1, i2=i2, row1=row1, row2=row2))
		} else {
			return(c())
		}
	}
	
	match3 <- function(i3, i1, i2, row1, row2) {
		val <- d3[i3,i1]
		if(is.na(val)) {
			return(c())
		} else {
			return(c(i1, i2, i3, row1[i2], row2[i3], val))
		}
	}
	
	results <- lapply(dimnames(d1)[[1]], match1)
	results <- unlist(results)
	if(!is.null(results)) {
		results <- as.data.frame(matrix(results, ncol=6, byrow=TRUE), stringsAsFactors=FALSE)
		for(i in 4:6) { results[,i] <- as.numeric(results[,i]) }
		invisible(results)
	} else {
		invisible(c())
	}
}
