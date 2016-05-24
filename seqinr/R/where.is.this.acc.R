where.is.this.acc <- function(acc, stopAtFirst = TRUE, ...){
  #
  # Argument check:
  #
  if(!is.character(acc)) stop("string expected for argument acc")
  #
  result <- character(0)
  #
  cat("Looking for available databases\n")  
  banks <- choosebank(...)
  nbanks <- length(banks)
  cat(paste("Looking for sequence with accession number", acc,
  "in the following ACNUC databases:\n"))
  print(banks)
  #
  # Looping over banks:
  #
  for(i in seq_len(nbanks)){
  	  cat(paste("\nTrying to open bank with name --->", banks[i], 
  	  "<--- ...", sep = ""))
    bkopenres <- try(choosebank(banks[i]))
    if(inherits(bkopenres, "try-error")){
      cat("... opening not OK, skipping this bank.\n")
  } else {
    cat("... and opening was OK.\n")
    cat(paste("==> Trying to find sequence", acc, "in bank", 
      banks[i], "..."))
    resquery <- try(query(".tmpquery", paste("AC=", acc)), silent = TRUE)
    if(inherits(resquery, "try-error")){
    	 cat("... not found here.\n")
    } else {
      cat("... *** FOUND *** here.\n")
      result <- c(result, banks[i])
    }
    closebank()
    if(length(result) != 0 && stopAtFirst){
      return(invisible(result))
    }
  }
}
#
# Print result summary:
#
  cat("\n\n")
  if(length(result) == 0){
  	  cat(paste("Sequence with accesion number", acc, 
  	  "was not found in available databases.\n Are you sure this is an accession number and not a sequence name?"))
  	} else {
    cat(paste("Sequence with accesion number", acc, 
    "was found in the following database(s):\n"))
    print(result)
  }
  invisible(result)
}

