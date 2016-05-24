Jaccard.Index <- function(x, y){
  xx <- as.vector(x)
  yy <- as.vector(y)

  if(length(xx) != length(yy)){
    return("Error: x and y have different numbers of elements")
  }

  id <- (!is.na(xx)) & (!is.na(yy))
  xx <- xx[id]
  yy <- yy[id]
  # xx[xx != 1] <- 0
  # yy[yy != 1] <- 0
  sum(xx == 1 & yy == 1) / sum(xx == 1 | yy == 1)
} # End of Jaccard.Index().
