strip.infraspec <- function(x, force.underscore = FALSE){
  if ( is.factor(x) ) x <- levels(x)[x]
  x <- gsub("_x_", "_x-", x) # handle times symbol in hybrids
  sepchar <- ifelse(length(grep("_", x)) != 0, "_", " ")
  x <- strsplit(x, sepchar)
  if ( force.underscore ) sepchar <- "_"
  x <- sapply(x, function(x, sc) paste(x[1:2], collapse = sc), sc = sepchar)
  x <- gsub("_x-", "_x_", x) # handle times symbol in hybrids
  x
}