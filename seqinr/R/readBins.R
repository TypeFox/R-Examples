readBins <- function(file, 
  colnames = c("allele.name", "size.bp", "minus.bp", "plus.bp")){
  src <- readLines(file)
  iPanel <- which(substr(src, start = 1, stop = 11) == "Panel Name\t")
  mycon <- textConnection(src[1:3])
  infos <- read.table(mycon, sep = "\t", fill = TRUE, header = FALSE)
  close(mycon)
  result <- list(infos = infos)

  starts <- iPanel + 1
  stops  <- c(iPanel[-1] - 1, length(src))
  for(i in seq_len(length(iPanel))){
  	  locsrc <- src[starts[i]:stops[i]]
  	  iMark <- which(substr(locsrc, start = 1, stop = 12) == "Marker Name\t")
  	  locres <- vector(mode = "list", length = length(iMark))
  	  locstarts <- iMark + 1
  	  locstops <- c(iMark[-1] - 1, length(locsrc))
  	  for(j in seq_len(length(iMark))){
  	  	  mycon <- textConnection(locsrc[locstarts[j]:locstops[j]])
  	  	  locres[[j]] <- read.table(mycon, sep = "", fill = TRUE) # changed
  	  	  colnames(locres[[j]])[1:length(colnames)] <- colnames
  	  	  close(mycon)
  	  	  names(locres)[j] <- unlist(strsplit(locsrc[iMark[j]], split = "\t"))[2]
  	  	}
  	  	#
  	  	# Check that the number of columns is 4 for all loci:
  	  	#
  	  	ncols <- sapply(locres, ncol)
  	  	if(any((ncols != 4))){
  	  		warning("A problem may have occur during importation")
  	  	}
    result[[i+1]] <- locres
    names(result)[i+1] <- unlist(strsplit(src[iPanel[i]], split = "\t"))[2]
  }
  return(result)
}
