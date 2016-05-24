parse.deconstruct <-
function(P) {
  i <- 1
  k <- 0
  if (length(P$sumset) == 0) {
    while (i <= length(P$children) & length(P$children) > 0 & length(P$divisor$children) > 0) {
      is.element <- FALSE
      for (j in 1:length(P$divisor$children)) {
        if (identical(P$children[[i]], P$divisor$children[[j]], attrib.as.set = TRUE)) {
          is.element <- TRUE
          k <- j
          break
        }
      }
      if (is.element) {
        P$divisor$children[[k]] <- NULL
        P$children[[i]] <- NULL
        i <- 0
      }
      i <- i + 1
    }
  }
  if (length(P$divisor$children) == 0) {
    P$fraction <- FALSE
    P$divisor <- list()
  }
  return(P)
}
