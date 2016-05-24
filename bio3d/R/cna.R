`cna` <-
function(cij, ...) {
  if(inherits(cij, "matrix")) {
    class(cij) <- c("matrix", "dccm")
    UseMethod("cna", cij)
  } else {
   class(cij) <- c("ensmb")
   UseMethod("cna", cij)
  }
}
