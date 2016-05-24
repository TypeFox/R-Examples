arrayRtoJSON <- function(
  Rarray
  ) {
  str = paste(Rarray, collapse="\",\"")
  str = paste("[\"", str, "\"]", sep="")
  str
}
