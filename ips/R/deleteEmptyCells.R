deleteEmptyCells <- function(DNAbin, 
                               nset = c("-", "n", "?"), 
                               quiet = FALSE){
  
  if ( !inherits(DNAbin, "DNAbin") ) 
    stop("'DNAbin' is not of class 'DNAbin'")
  
  ## convert character to raw
  nn <- as.raw(c(240, 2, 4))
  names(nn) <- c("n", "?", "-")
  nset <- nn[nset]
	
  ## function that detects non-empty strings
	isNotEmpty <- function(x, nset){
		ifelse(all(unique(x) %in% nset), FALSE, TRUE)
	}
  
  size <- dim(DNAbin)
	
	## columns
	colind <- which(apply(DNAbin, 2, isNotEmpty, nset = nset))
	DNAbin <- DNAbin[, colind]
	
	## rows
	rowind <- which(apply(DNAbin, 1, isNotEmpty, nset = nset))
	DNAbin <- DNAbin[rowind, ]
  
  ## screen output (if desired)
  if ( !quiet ) {
    size <- size - dim(DNAbin)
    cat("\n\t", size[1], "rows deleted from alignment")
    cat("\n\t", size[2], "columns deleted from alignment")
  }  
	DNAbin
}

