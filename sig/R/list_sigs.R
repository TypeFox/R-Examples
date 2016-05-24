#' List the signatures of all functions
#' 
#' Lists the signatures of all functions in an environment
#' or file.
#' @param x An environment or the the path to a file.
#' @param pattern An optional regular expression. Only names
#' matching pattern are returned.       
#' @param ... Currently ignored
#' @return An object of class \code{siglist}, which is a list
#' of \code{sig} obejcts.  
#' @examples
#' #From a package
#' list_sigs(pkg2env(graphics))
#' #Just functions beginning with 'a'.
#' list_sigs(pkg2env(graphics), pattern = "^a")
#' #From a file
#' list_sigs(system.file("extdata", "sample.R", package = "sig"))
#' @export
list_sigs <- function(x, pattern = NULL, ...)
{
  UseMethod("list_sigs")  
}
  
#' @rdname list_sigs  
#' @method list_sigs default
#' @export
list_sigs.default <- function(x, ...)
{
  x <- as.environment(x)
  NextMethod("list_sigs")
}

#' @rdname list_sigs  
#' @method list_sigs environment
#' @export 
list_sigs.environment <- function(x, pattern = NULL, ...)
{
  fns <- Filter(is.function, as.list(x))
  if(!is.null(pattern))
  {
    fns <- fns[grepl(pattern, names(fns))]
  }
  o <- order(names(fns))
  fns <- fns[o]
  structure(
    mapply(
      sig,
      fns,
      fix_fn_names(names(fns)),
      SIMPLIFY = FALSE
    ),
    class = c("siglist", "list")
  )
}

#' @rdname list_sigs  
#' @method list_sigs character
#' @export       
list_sigs.character <- function(x, ...)
{
  e <- source_to_new_env(x)
  list_sigs(e, ...)
}

#' @rdname print.siglist
#' @method toString siglist
#' @export
toString.siglist <- function(x, width = getOption("width"), ...)
{
  unlist(
    lapply(
      x, 
      function(s) c(toString(s, width = width, ...), "")
    ),
    use.names = FALSE
  )  
}

#' Print a siglist object
#' 
#' Prints a list of function signature objects.
#' 
#' @method print siglist
#' @param x An object of class \code{siglist}.
#' @param width Width of string to display.
#' @param ... Passed to the equivalent \code{sig} method.
#' @return \code{toString} creates a string representation of a function signature. 
#' \code{print} is mostly invoked for the side effect of printing a function 
#' signature, invisibly returning its input.
#' @examples
#' method_sigs <- list_sigs(pkg2env(methods))
#' print(method_sigs)
#' print(method_sigs, width = 40)
#' print(method_sigs, width = 40, exdent = 2)
#' toString(method_sigs)
#' @export
print.siglist <- function(x, width = getOption("width"), ...)
{
   print_engine(x, width,  ...)
}
