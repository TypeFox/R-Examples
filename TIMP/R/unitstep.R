"unitstep" <- function(x) {
  xx <- x
  xx[which(x>=0)] <- 1
  xx[which(x<0)] <- 0 
  xx
}
