#' Private class
#' 
#' This is a virtual class to be contained in other class definitions. It overrides the default subset functions \code{$} and \code{[[} such that private member of a class can not be accessed. Private is every object which has a name with a leading "." (\code{grepl("^\\\\.", name)}). After this check the standard method for class 'envRefClass' is called or an error is reported.
#' 
#' @seealso \link{defineRefClass}
#' @export
#' @rdname Private
#' @examples
#' ClassWithPrivateField <- defineRefClass({
#'   Class <- "ClassWithPrivateField"
#'   contains <- "Private"
#'   
#'   .p <- "numeric"
#'   
#'   getP <- function() .p
#'   setP <- function(v) .self$.p <- v
#' })
#' 
#' test <- ClassWithPrivateField()
#' stopifnot(inherits(try(test$.p), "try-error"))
#' stopifnot(inherits(try(test$.p <- 2), "try-error"))
#' stopifnot(inherits(try(test[[".p"]]), "try-error"))
#' stopifnot(inherits(try(test[[".p"]] <- 2), "try-error"))
setRefClass("Private", contains = "VIRTUAL")

#' @export
#' @rdname Private
#' 
#' @param x the object
#' @param name name of field or method
setMethod("$", "Private", function(x, name) {
  
  .self <- selectMethod("$", "envRefClass")(x, ".self")
  
  callFromInside <- any(sapply(envirSearch(list(parent.frame())), 
                               identical, y = as.environment(.self)))

  if(!callFromInside & grepl("^\\.", name)) stop("Restricted access!")
  
  selectMethod("$", "envRefClass")(x, substitute(name))
})

#' @export
#' @rdname Private
#' 
#' @param value any object
setMethod("$<-", "Private", function(x, name, value) {
  
  .self <- selectMethod("$", "envRefClass")(x, ".self")
  
  callFromInside <- any(sapply(envirSearch(list(parent.frame())), 
                               identical, y = as.environment(.self)))
  
  if(!callFromInside & grepl("^\\.", name)) stop("Restricted access!")
  
  selectMethod("$<-", "envRefClass")(x, substitute(name), value)
})

#' @export
#' @rdname Private
#' 
#' @param i like name
#' @param j ignored
#' @param ... ignored
setMethod("[[", "Private", function (x, i, j, ...) {
  stop("This method is disabled for this class!")
})

#' @export
#' @rdname Private
setMethod("[[<-", "Private", function (x, i, j, ..., value) {
  stop("This method is disabled for this class!")
})
