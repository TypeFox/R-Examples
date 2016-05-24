formatPercent <- function(x, digits = 1){
  res <- paste(formatC(x * 100, digits, format = "f"), "%", sep = "")
  return(res)
}
