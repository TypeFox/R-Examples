trimEnds <- function(x, min.n.seq = 4){
  
  if ( !inherits(x, "DNAbin") ){
    stop("'x' is not of class 'DNAbin'")
  }
  if ( !is.matrix(x) ){
    stop("'x' must be a matrix")
  }
  
  ## internal function definition
  replaceWithN <- function(x){
    id <- x == as.raw(4)
    if ( length(id) > 0 & any(id[c(1, length(id))])){
      id <- which(id)
      getIndex <- function(x){
        
        ## indices of head to be filled with N
        for (i in seq_along(id) - 1) {
          if ( any(id[1:(i + 1)] != (1:(i + 1))) ) break
        }
    
        ## indices of tail to be filled with N
        id <- rev(id)
        jj <- head(id, 1); j <- jj - 1
        for (k in seq_along(id)[-1] ) {
          if ( any(id[1:k] != (jj:j)) ) break
          j <- j -1
        }
        j <- j + 1
        id <- c(0:i, j:jj)
        id[id != 0]
      }
      id <- getIndex(id)
      x[id] <- as.raw(240)
    } 
    return(x)
  }
  
  ## function body
  x <- t(apply(x, 1, replaceWithN))
  class(x) <- "DNAbin"
  
  ## trim ends to 'min.n.seq' bases
  ## ------------------------------
  b <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b){
    length(which(x %in% b))
  }
  m <- apply(x, 2, percentInformation, b)
  if ( max(m) < min.n.seq ) stop("alignment contains less sequences then required")
  m <- range(which(m >= min.n.seq))
  m <- seq(from = m[1], to = m[2])
  x <- x[, m]
  
  x
}



