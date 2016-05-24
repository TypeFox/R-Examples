"as.fasta" <- function(x, id=NULL, ...) {
  cl <- match.call()
  
  if(is.list(x)) {
    if(is.null(id))
      id <- x$id
    
    x <- x$ali
  }
  
  if(is.vector(x)) {
    if(any(nchar(x)>1))
      stop("provide a matrix/vector of one letter amino acid codes")
      ##x <- seqbind(lapply(lapply(x, strsplit, ""), unlist))

    x <- as.matrix(t(x))
  }
  
  if(is.matrix(x)) {
    if(is.null(id))
      id <- rownames(x)
    
    if(is.null(id))
      id <- paste("seq",1:nrow(x), sep="")
    
    if(any(id=="") | any(is.na(id))) {
      id[id==""] <- NA
      inds <- which(is.na(id))
      id[inds] <- paste("seq", inds, sep="")
    }
    
    if(nrow(x) != length(id))
      stop("length of 'id' does not match number of rows in alignment")
    
    rownames(x) <- id
  }

  else {
    stop("provide a sequence character matrix/vector")
  }
  
  out <- list(id=id, ali=x, call=cl)
  class(out) <- "fasta"
  return(out)
  
}
