##
## Utility functions for handling examples in book
##
## List names of examples
listEx <- function(){
  p <- file.path(find.package("FRAPO"), "BookEx")
  f <- list.files(p, pattern = "*.R")
  ans <- sub("^([^.]*).*", "\\1", f)
  ans
}
## Display example
showEx <- function(Example){
  Exs <- listEx()
  if(!(all(Example %in% Exs))){
    stop("\nName of example is not included.\n")
  }
  Ex <- paste(Example, "R", sep = ".")
  f <- file.path(find.package("FRAPO"), "BookEx", Ex)
  ans <- file.show(f)
  invisible(ans)
}
## Save a copy of a example in the working directory
saveEx <- function(Example){
  Exs <- listEx()
  if(!(all(Example %in% Exs))){
    stop("\nName of example is not included.\n")
  }
  Ex <- paste(Example, "R", sep = ".")
  f <- file.path(find.package("FRAPO"), "BookEx", Ex)
  file.copy(f, getwd(), copy.mode = TRUE)
}
## Edit a copy of an example in the working directory
editEx <- function(Example, ...){
  Ex <- paste(Example, "R", sep = ".")
  saveEx(Example)
  file.edit(Ex, ...)  
}
## Run an example
runEx <- function(Example, ...){
  Exs <- listEx()
  if(!(all(Example %in% Exs))){
    stop("\nName of example is not included.\n")
  }
  Ex <- paste(Example, "R", sep = ".")
  f <- file.path(find.package("FRAPO"), "BookEx", Ex)
  source(f, ...)
}
