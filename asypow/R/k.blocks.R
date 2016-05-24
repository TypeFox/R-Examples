k.blocks <- function( ... ) {

      blocks <- list(...)

      n <- 0
      for (block in blocks) n <- n + dim(block)[1]
      bmat <- matrix(0, n, n)

      ind <- 1
      for (block in blocks) {
           endind <- ind + dim(block)[1] - 1
           bmat[ind:endind,ind:endind] <- block
           ind <- endind + 1
      }

      return(bmat)
}
