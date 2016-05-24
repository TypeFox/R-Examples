### To unlist a strsplit object.
unstrsplit <- function(x, y, ...){
  unlist(strsplit(x, y, ...))
} # End of unstrsplit().


### Functions for summary() and print.summary().
my.format <- function(x, digits = max(4, getOption("digits") - 3)){
  paste(formatC(x, format = "f", width = -1, digits = digits), collapse = " ")
} # End of my.format().

my.cat <- function(...){
  cat(..., sep = "")
} # End of my.cat().

my.print <- function(x, digits = max(4, getOption("digits") - 3)){
  print.default(x, na.print = "", quote = FALSE, right = TRUE, digits = digits)
} # End of my.print().

