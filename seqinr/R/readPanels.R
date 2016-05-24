readPanels <- function(file,
  colnames = c("marker", "dye.col", "min.bp", "max.bp", "exp.pcg", "repeat.bp",
    "stutter.pc", "uknw", "allele names")){
  src <- readLines(file)
  iPanel <- which(substr(src, start = 1, stop = 6) == "Panel\t")

  infos <- src[1:(iPanel[1] - 1)]
  result <- list(infos = infos)

  starts <- iPanel + 1
  stops  <- c(iPanel[-1] - 1, length(src))
  
  for(i in seq_len(length(iPanel))){
  	 toimport <- src[starts[i]:stops[i]]
  	 # Change all runs of tabulations by a single one:
  	 toimport <- gsub("\t{2,}", "\t", toimport)
	 mycon <- textConnection(toimport)
    result[[i+1]] <- read.table(mycon, sep = "\t", quote = "")
    close(mycon)
    # remove empty columns
    # (I guess this is no more necessary now that runs of
    # tabulation are preprocessed)
    tokeep <- rep(TRUE, ncol(result[[i+1]]))
    for(j in 1:ncol(result[[i+1]])){
      if(all(is.na(result[[i+1]][,j]))){
      	tokeep[j] <- FALSE
      }
    }
    result[[i+1]] <- result[[i+1]][, tokeep]
    # There must be 9 columns
    if(ncol(result[[i+1]]) != 9) stop("wrong column number")
    colnames(result[[i+1]]) <- colnames
    
    headeritems <- unlist(strsplit(src[iPanel[i]], split = "\t"))
    # remove empty elements
    if(any(nchar(headeritems) == 0)){
      headeritems <- headeritems[-which(nchar(headeritems) == 0)]
    }
    names(result)[i+1] <- headeritems[2]
  }
  return(result)
}
