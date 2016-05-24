# Nov 2014: add option to retain dimensionality on request
# except,... for arbitrary rank array, what would the user want? 
# limit the fancy alternative layout to 2D only 
#
short <- function (x = seq(1, 20), numel = 4, skipel = 0, ynam = deparse(substitute(x)), dorows=FALSE)  {
# check for array
	if(is.array(x) && length(dim(x)) >2 && dorows) {
		warning("Cannot return layers of an array. Treating as vector")
		dorows <- FALSE 
	}
	if(is.null(dim(x)) &&dorows) warning("input is not array; dorows=TRUE is ignored.")
    ynam <- as.character(ynam)
    ynam <- gsub(" ", "", ynam)
    if (is.list(x)) 
        x <- unlist(t(x))
    if (2 * numel >= length(x)) {
        print(x)
    }
    else {
# new code for matrices
		if(!is.null(dim(x))&& dorows && length(dim(x)==2) ) {
			# numel is now numrow; but check that it's not too big
			if(2*numel>=nrow(x)) {
				print(x)
			}else {
#   check skipel - if sum of skipel and numel exceeds nrow, then just
# return the last numel as the"top" and the first numel rows as the "bottom"
				skipel <- min(skipel,nrow(x)-numel)
				frist = 1 + skipel
	        	last = numel + skipel
			# make row names and that's all we need
				rownames(x)<-1:nrow(x)
				print(x[c( frist:last,  ((nrow(x)-skipel-numel+1):(nrow(x)-skipel))) ,] ) 
			}
		} else { 
#in this case, don't have a matrix, or dorows==FALSEsh
# here, clamp last to length(vectorized input) and frist to numel less than that.	
			skipel <- min(skipel,length(x)-numel)
	        frist = 1 + skipel
	        last = numel + skipel
	        cat(paste(ynam, "[", frist, "] thru ", ynam, "[", last, 
	            "]\n", sep = ""))
	        print(x[frist:last])
	        cat(" ... \n")
	 
	        cat(paste(ynam, "[", length(x) - numel - skipel + 1, 
	            "] thru ", ynam, "[", length(x) - skipel, "]\n", 
	            sep = ""))
	        print(x[(length(x) - numel - skipel + 1):(length(x) - 
	            skipel)])
		}
    }
}
