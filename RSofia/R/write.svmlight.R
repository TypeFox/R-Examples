write.svmlight <- function(labels, data, file, ...) {

  if(!is.vector(labels))
    stop("labels must be a vector")

  if(!is.matrix(data))
    stop("data must be a matrix")

  val <- .Call("svmlight_writer", file, data, labels, PACKAGE = "RSofia")

}
