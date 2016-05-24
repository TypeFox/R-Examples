### Generate a (column-wised) list structure for given # of rows and columns.
my.generate.list <- function(value, nrow, ncol){
  ret <- vector(mode = "list", length = ncol)
  for(i in 1:ncol){
    ret[[i]] <- rep(value, length = nrow)
  }
  ret
} # End of my.generate.list()

