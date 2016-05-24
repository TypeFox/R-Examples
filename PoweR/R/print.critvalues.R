print.critvalues <- function(x, digits = 3, latex.output = FALSE, template = 1, ...) {

  class(x) <- paste("critvalues",template,sep="")

  print(x, digits, latex.output, ...)

}
