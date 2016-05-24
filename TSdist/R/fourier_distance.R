

# This function calculates a distance based on the Fourier Coefficients
FourierDistance <- function(x, y, n = (floor(length(x) / 2) + 1)) {

  err <- try(FInitialCheck(x, y, n))
  if (class(err) == "try-error") {return(NA)}
    
  # The discrete fourier transform is calculated for both series.
  fft1 <- fft(x)
  fft2 <- fft(y)
    
  # The distance is calculated by using the first n coefficients and applying the 
  # Euclidean distance between them.
  d <- sqrt(sum(Mod(fft1[1:n] - fft2[1:n]) ^ 2))

  return(d)
}



# This function checks for possible initial errors: 
FInitialCheck <- function(x, y, n) { 
  
  if (! is.numeric(x) | ! is.numeric(y)) {
    stop('The series must be numeric', call.=FALSE)
  }
  if (! is.vector(x) | ! is.vector(y)) {
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (length(x) <= 1 | length(y) <= 1) {
    stop('The series must have a more than one point', call.=FALSE)
  }
  if (length(x) != length(y)) {
    stop('The length of the series must be equal', call.=FALSE)
  }
  if (n > length(x)) {
    stop('The number of Fourier Coefficients included may not exceed 
         the length of the series', call.=FALSE)
  }
  if (n > (length(x) / 2 + 1)) {
    warning('The number of coefficients considered is larger than necessary.', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))) {
    stop('There are missing values in the series', call.=FALSE)
  } 
}
