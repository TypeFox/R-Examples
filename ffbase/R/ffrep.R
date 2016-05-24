#' Replicate Elements of \code{ff} vectors.
#'
#' Similar as \code{rep.int} in the base package but for \code{ff} vectors.
#'
#' @export
#' @example ../examples/ffrep.R
#' @param x an integer \code{ff} vector
#' @param times integer \code{ff} vector giving the (non-negative) number of times to repeat each element if of length length(x), 
#' or an integer of length 1 indicating how many times to to repeat the whole vector. Negative or NA values are an error.
#' @return An ff vector of integers where x is recycled 
#' @seealso \code{\link[base]{rep.int}}
ffrep.int <- function(x, times){		
	if(length(times) == 1){
		result <- ff(x, length=length(x)*times[])
	}else{
		stopifnot(length(times) == length(x))
		BATCHBYTES <- getOption("ffbatchbytes")
		## Compute how to do the batched replication (in memory is times x .rambytes = the .rambytes of the times itself
		rambytes <- times * .rambytes[vmode(x)] + .rambytes[vmode(times)]
		## find indexes of the locations of the chunks
		idxlocations <- integer(0)
		i.last <- 0
	 	for (i in chunk(rambytes)){
 			overflows <- grouprunningcumsumindex(x=rambytes[i], max=BATCHBYTES, currentcumul=i.last)
 			idxlocations <- append(idxlocations, overflows$overflowidx + min(i))
 			i.last <- overflows$currentcumul
 		}	
 		idxlocations <- idxlocations-1
 		idxlocations <- unique(c(1, idxlocations, length(rambytes)))
 		## make a list of ri to do the batchwise rep
 		n <- length(idxlocations)
 		s <- seq_len(n-1)
 		ret <- list()
 		for (i in s) {
 			if(i == 1){ 	
 				ret[[i]] <- ri(idxlocations[i], idxlocations[i+1])			
 			}else{
 				ret[[i]] <- ri(idxlocations[i]+1, idxlocations[i+1])
 			}    	
    }
    ##
    ## do the replication
    ##
		result <- ff::clone.ff(x, initdata=NULL, length=sum(times, na.rm=TRUE))	
		from <- 1
		for (i in ret){
			## extract the ff data and replicate
    	toupdate <- rep.int(x[i], times[i])    	
    	to <- from+length(toupdate)-1
    	iri <- ri(from = from, to = to)
    	## update the result ff
    	result[iri] <- toupdate
    	from <- to+1
  	}
	}
	result
}

