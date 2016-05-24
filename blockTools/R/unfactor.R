unfactor <- function(x){
  tmp <- x
  tmp2 <- as.numeric(levels(tmp)[as.integer(tmp)])
  return(tmp2)
}
