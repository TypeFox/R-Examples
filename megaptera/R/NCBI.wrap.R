NCBI.wrap <- function(x, BOOLEAN = "OR", field = "all"){
  x <- paste(x, "[xxx]", sep = "")
  x <- paste(x, collapse = paste(rep(" ", 2), collapse = BOOLEAN))
  x <- paste("(", x, ")", sep = "")
  xx <- vector()
  for (i in field){
    xx <- c(xx, gsub("xxx", i, x))
  }
  xx
}