is.temporal <- function(x) {
  return(length(grep("^[a-z]?[0-9]{4}$",x))==length(x))
}