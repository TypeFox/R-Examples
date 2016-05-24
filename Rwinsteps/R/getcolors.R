getcolors <- function(col, len) {

  if(length(col) == 1) {
    if(tolower(col) == "r")
      col <- rainbow(len)
  }
  return(rep(col, length.out = len))
}
