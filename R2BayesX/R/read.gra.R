read.gra <- function(file, sorted = FALSE, sep = " ")
{    
  ## read information from the graph file
  gfile <- strsplit(readLines(file), sep)

  ## get region names and neighbors
  if(length(gfile[[2]]) < 2) {
    ids <- rep(1:3, length = length(gfile) - 1)
    regions <- unlist(gfile[c(2:length(gfile))[ids == 1]])
    neighbors <- gfile[c(2:length(gfile))[ids == 3]]
    foo <- function(x) as.integer(x) + 1
  } else {
    regions <- sapply(gfile[-1], function(x) x[1])
    neighbors <- sapply(gfile[-1], function(x) x[-c(1:2)])
    foo <- function(x) as.integer(x)
  }
  neighbors <- lapply(neighbors, foo)

  ## extract the number of regions
  n <- length(regions)
  cat("Note: map consists of", n, "regions\n")

  cat("Creating adjacency matrix ...\n")

  ## create the adjacency matrix
  adjmat <- matrix(0, n, n)
  
  ## number of neighbors for each region
  diag(adjmat) <- sapply(neighbors, length)

  ## loop over the regions
  for(i in 1:n) {
    ## off-diagonal entries in adjmat
    ## Note: the indices are zero-based!
    if(length(neighbors[[i]]))
      adjmat[i, neighbors[[i]]] <- -1
  }

  colnames(adjmat) <- rownames(adjmat) <- regions

  cat("finished\n")
    
  if(sorted) {
    ## sort regions and adjacency matrix
    if(sum(is.na(as.numeric(regions))) == 0){
      regions <- as.numeric(regions)
      cat("Note: regions sorted by number\n") 
    } else cat("Note: regions sorted by name\n")
    ord <- order(regions)
    regions <- sort(regions)
    adjmat <- adjmat[ord, ord]
  }

  rownames(adjmat) <- colnames(adjmat) <- regions

  class(adjmat) <- c("gra", "matrix")
  return(adjmat)
}

