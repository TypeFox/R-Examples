StripTrailingSlashes <- function(x) {
  while( grepl("[/\\\\]$", x) ) x <- substring(x, 1, nchar(x) - 1)
  return(x)
}