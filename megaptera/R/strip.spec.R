strip.spec <- function(x){
  if ( is.factor(x) ) x <- levels(x)[x]
  sepchar <- ifelse(length(grep("_", x)) != 0, "_", " ")
  x <- strsplit(x, sepchar)
  sapply(x, function(x) x[1])
}