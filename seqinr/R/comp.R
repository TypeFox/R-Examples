##########################################
# complements a sequence
###########################################

comp <- function(seq, forceToLower = TRUE, ambiguous = FALSE){
  
  if(all(seq %in% LETTERS)){
  	  isUpper <- TRUE
  	} else {
  	  isUpper <- FALSE
  	}
  	
  	seq <- tolower(seq)

  result <- as.vector(n2s((3-s2n(seq))))
  #
  # More work is required if ambiguous bases are taken into account
  #
  
  if(ambiguous){
    result[which(seq == "b")] <- "v"
    result[which(seq == "d")] <- "h"
    result[which(seq == "h")] <- "d"
    result[which(seq == "k")] <- "m"
    result[which(seq == "m")] <- "k"
    result[which(seq == "s")] <- "s"
    result[which(seq == "v")] <- "b"
    result[which(seq == "w")] <- "w"
    result[which(seq == "n")] <- "n"    
    result[which(seq == "y")] <- "y"
    result[which(seq == "r")] <- "r"
  }
  
  # Checking for N in the sequence, thanks to Jeremy Shearman. 
  
  result[which(seq == "n")] <- "n"	
  
  
  if(isUpper && !forceToLower){
  	  result <- toupper(result)
  	}
  return(result)
}

