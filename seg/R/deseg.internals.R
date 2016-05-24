# ------------------------------------------------------------------------------
# Internal functions used by 'deseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
.decomp <- function(data) {

  dx <- apply(data, 1, sum) # total population in each row (i.e., area)
  removeID <- which(dx == 0)
  removeL <- length(removeID)
  if (removeL > 0) {
    warning("remove ", removeL, " rows with no population", call. = FALSE)
    dx <- dx[-removeID]
    data <- data[-removeID,]
  }  

  sx <- data / dx
  sx <- apply(sx * sx, 1, sum)
  sum(dx * sx) / sum(dx)
}

.decompL <- function(data) {
  
  groupsize <- apply(data, 2, sum)
  numpoints <- nrow(data)
  
  tmp <- rep(groupsize / numpoints, numpoints)
  dataL <- matrix(tmp, nrow = numpoints, byrow = TRUE)
  
  .decomp(dataL)
}

.decompC <- function(data) {
  
  tmp <- sum(data) / (nrow(data) * ncol(data))
  dataC <- matrix(tmp, nrow = nrow(data), ncol = ncol(data))
  
  .decomp(dataC)
}
