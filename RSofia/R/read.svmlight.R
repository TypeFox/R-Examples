read.svmlight <- function(file) {

  x <- .Call("svmlight_reader", file, PACKAGE = "RSofia")
  
  return (x)

}
