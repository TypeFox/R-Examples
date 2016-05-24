k.blocks.info <- function( blocks ) {

      n <- dim(blocks[[1]])[1]
      size <- length(blocks) * n
      info <- matrix(0, size, size)

      ind <- 1
      for (block in blocks) {
           endind <- ind + n - 1
           info[ind:endind,ind:endind] <- block
           ind <- endind + 1
      }

      return(info)
}
