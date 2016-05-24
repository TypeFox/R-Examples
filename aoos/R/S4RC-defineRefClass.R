#' Define a Reference Class
#'
#' This is a wrapper around \code{\link{setRefClass}}. All arguments are defined in an expression (instead of lists) which improves readability of the code. Besides that, no additional features are added.
#' 
#' @param expr an expression
#'
#' @seealso \link{Private-class}
#' @export
#' @examples
#' \dontrun{
#'   vignette("Introduction", "aoos")
#' }
#'
#' # Minimal example:
#' Test <- defineRefClass({
#'   Class <- "Test" # this is passed as argument to setRefClass
#'   x <- "character" # all objects which are not functions are fields
#'   do <- function() cat("Yes, Yes, I'm working...") # a method
#' })
#'
#' test <- Test()
#' test$x <- "a"
#' test$do()
#'
#' # Inheritance and privacy:
#' pTest <- defineRefClass({
#'   Class <- "pTest"
#'   # Privacy is solved by inheriting from a class 'Private' which redefines
#'   # the methods for access.
#'   contains <- c("Test", "Private") # passed as argument to setRefClass
#'
#'   .y <- "numeric" # this is going to be 'private'
#'
#'   doSomething <- function() {
#'     .self$.y <- 42
#'     cat(x, .y, "\n")
#'     invisible(.self)
#'   }
#' })
#'
#' instance <- pTest()
#' instance$x <- "Value of .y:"
#' instance$doSomething()
#'
#' # A notion of privacy:
#' stopifnot(inherits(try(instance$.y), "try-error"))
#' stopifnot(inherits(try(instance$.y <- 2), "try-error"))
defineRefClass <- function(expr) {

  mc <- match.call()
  e <- new.env()
  eval(mc$expr, e)
  argsList <- as.list(e, all.names = TRUE)

  contains <- as.character(argsList$contains)
  argsList$contains <- NULL
  Class <- as.character(argsList$Class)
  argsList$Class <- NULL

  argsList <- combineListElements(
    argsList,
    sapply(argsList, Negate(function(e) inherits(e, what = "function"))),
    "fields")

  argsList <- combineListElements(
    argsList,
    sapply(argsList, inherits, what = "function"),
    "methods")

  argsList$contains <- contains
  argsList$Class <- Class
  argsList$where <- parent.frame()

  do.call(setRefClass, args = argsList)

}

combineListElements <- function(l, ind, name) {
  ind <- as.logical(ind)
  newElement <- l[ind]
  l[ind] <- NULL
  l[[name]] <- newElement
  l
}
