is.spatial <- function(x) {
  return(length(grep("^(([A-Z]{3})|(glob)|([A-Z]+[\\._][0-9]+))$",x))==length(x))
}



