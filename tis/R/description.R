description <- function(x)
  attr(x, "description")

"description<-" <- function(x, value){
  if(is.character(value) && nchar(value) == 0) value <- NULL
  attr(x, "description") <- value
  invisible(x)
}

documentation <- function(x)
  attr(x, "documentation")

"documentation<-" <- function(x, value){
  if(is.character(value) && nchar(value) == 0) value <- NULL
  attr(x, "documentation") <- value
  invisible(x)
}
