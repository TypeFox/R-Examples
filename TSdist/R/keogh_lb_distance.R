LBKeoghDistance <- function(x, y, window.size) {
  
  if (class(try(lbKeoghInitialCheck(x, y, window.size))) == "try-error") {
    return(NA)
  } else {
  
  # The upper envelope of x is built.
  upper.env <- rollapply(x, window.size, max, partial=TRUE)
  # The lower envelope of y is built.
  lower.env <- rollapply(x, window.size, min, partial=TRUE)
  # The distance is calculated:
  D <- c(1:length(y)) * 0
  ind1 <- which(y > upper.env)
  D[ind1] <- (y[ind1] - upper.env[ind1]) ^ 2
  ind2 <- which(y < lower.env)
  D[ind2] <- (y[ind2] - lower.env[ind2]) ^ 2
  d <- sqrt(sum(D))
  return(d)
  }
}



# This function checks for possible initial errors: 
lbKeoghInitialCheck <- function(x, y, window.size) {
  
  if (! is.numeric(x) | ! is.numeric(y)) {
    stop('The series must be numeric', call.=FALSE)
  }
  if (! is.vector(x) | ! is.vector(y)) {
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (length(x) < 1 | length(y) < 1) {
    stop('The series must have at least one point', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))) {
    stop('There are missing values in the series', call.=FALSE)
  } 
  if (length(x) != length(y)) {
    stop('The series must have the same length', call.=FALSE)
  }
  if (window.size %% 2 != 1) {
    stop('For the Sakoe-Chiba band, the window must be symmetric 
         and window.size must be even', call.=FALSE)
  } 
  if (window.size > length(x)) {
    stop('The width of the window should not exceed the length 
         of the series', call.=FALSE)
  } 
  
}