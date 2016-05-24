# ------------------------------------------------------------------------------
# Function 'seg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# Last update: 22 December 2013
# Depends: -
# ------------------------------------------------------------------------------
seg <- function(data, nb) {

  message("Note: seg() will be deprecated soon.")
  message("Use dissim() instead. See help(dissim) for details.")
  
  if (ncol(data) > 2) {
    warning("'data' has more than two columns; only the first two are used",
            call. = FALSE)
    data <- data[,1:2]
  }
  if (any(data < 0))
    stop("negative value(s) in 'data'", call. = FALSE)
  colsum <- apply(data, 2, sum)
  if (any(colsum <= 0))
    stop("the sum of each column in 'data' must be > 0", call. = FALSE)
      
  # Duncan and Duncan's index of dissimilarity
  b <- data[,1]/sum(data[,1])     # Blacks
  w <- data[,2]/sum(data[,2])     # Whites
  d <- sum(abs(b-w))/2

  if (!missing(nb)) {     # If information on geographic distance between
                          # spatial units is provided:
    if (!is.matrix(nb))
      stop("'nb' must be a matrix", call. = FALSE)
    else if (nrow(nb) != ncol(nb))
      stop("'nb' must be a square matrix", call. = FALSE)
    else if (nrow(nb) != nrow(data))
      stop("nrow(nb) must match nrow(data)", call. = FALSE)
    
    if (sum(nb) != 1)
      warning("the sum of all elements in 'nb' does not equal 1", call. = FALSE)

    rowsum <- apply(data, 1, sum)
    removeID <- which(rowsum == 0)
    removeL <- length(removeID)
    if (removeL > 0) {
      warning("remove ", removeL, " rows with no population", call. = FALSE)
      rowsum <- rowsum[-removeID]
      data <- data[-removeID,]
      nb <- nb[-removeID, -removeID]
    }
    
    # Black proportions in census tracts
    z <- data[,1] / rowsum
    # Additional spatial component value
    spstr <- 0
    nbvec <- as.vector(nb)
    INDEX <- which(nbvec != 0)
    for (i in 1:length(INDEX)) {
      rowID <- INDEX[i] %% nrow(nb)
      colID <- INDEX[i] %/% nrow(nb)
      if (rowID == 0)
        rowID <- nrow(nb)
      else
        colID <- colID + 1
      spstr <- spstr + (abs(z[colID] - z[rowID]) * nbvec[INDEX[i]])
    }
    d <- d - spstr
  }
  
  as.vector(d)
}
