kaks <- function(x, verbose=FALSE, debug = FALSE, forceUpperCase = TRUE){
    #
    # Check argument class:
    #
    if(attr(x,"class") != "alignment") stop("object x must be of class alignment")
    if(debug){
      cat("<--- Argument x storage is --->\n")
      print(str(x))
      cat("<--- Argument x storage is --->\n")
    }
    #
    # Check that there are at least two sequences in the alignment:
    #
    if(x$nb < 2){
      warning("there should be at least two sequences in the alignment")
      return(NA)
    }
    #
    # Check that all sequences are of the same length:
    #
    lseqs <- nchar(x$seq)
    if( !all(lseqs == lseqs[1])) {
      warning("all sequences should be the same length in an alignment")
      return(NA)
    }
    #
    # Check that the length of sequences is a mutiple of 3 since we are dealing
    # with coding sequences here:
    #
    if( lseqs[1] %% 3 != 0){
      warning("sequence lengths are not a multiple of 3")
      return(NA)
    }
    #
    # Force sequences characters to upper case letters when at least one
    # one 'a', 'c', 'g', or 't' is found in the sequences:
    #
    if(forceUpperCase){
      if( length(grep("[acgt]", x$seq)) != 0){
        x$seq <- toupper(x$seq)
      }
    }
    #
    # Call internal C function:
    #
    l <- .Call("kaks", x$seq, x$nb, debug, PACKAGE = "seqinr")
    if(debug){
      cat("<--- Result l storage is --->\n")
      print(str(l))
      print (l)
      cat("<--- Result l storage is --->\n")
    }
    #
    # If the sequences names are missing, we call them seq1, seq2, and so on:
    #
    if( is.null(x$nam) ) x$nam <- paste("seq", seq_len(x$nb), sep = "")
    
    #
    # This is to compute the list of results:
    #
    mkresult <- function(k){
      tmp <- matrix( k, x$nb, x$nb, byrow = TRUE, dimnames = list(x$nam, x$nam))
      as.dist(t(tmp))
    }
    #result <- lapply(l[seq_len(4)], mkresult)
    if (verbose)	
    	{
	result <- lapply(l[seq_len(13)], mkresult)
	check <- result[[5]] + result[[6]] +  result[[7]]
	result[[14]] <-check
	names(result) <- c("ka", "ks", "vka", "vks","l0","l2","l4","a0","a2","a4", "b0","b2","b4","checksuml")
	return(result)
	 }
    else
    	{
	result <- lapply(l[seq_len(4)], mkresult)
	names(result) <- c("ka", "ks", "vka", "vks")
    	return(result)
	}
}

