#' Summarise function complexity of a file or environment
#' 
#' @param x A path to an R file or an environment.
#' @param too_many_args Upper bound for a sensible number of args.
#' @param too_many_lines Upper bound for a sensible number of lines.
#' @param ... Passed to \code{sig_report.environment}.
#' @return An object of class ``sigreport'' with the elements.
#' \itemize{
#'   \item{n_vars}{Number of variables.}
#'   \item{n_fns}{Number of functions.}                            
#'   \item{n_args}{Table of the number of args of each function.}
#'   \item{too_many_args}{Upper bound for a sensible number of args.}
#'   \item{fns_with_many_args}{Names of each function with more args 
#'     than \code{too_many_args}.}
#'   \item{n_lines}{Table of the number of lines of each function body.} 
#'   \item{too_many_lines}{Upper bound for a sensible number of lines.}
#'   \item{long_fns}{Names of each function with more lines than 
#'   \code{too_many_lines}.}
#' }
#' @details
#' \code{sig_report} summarises the number of input arguments and the 
#' number of lines of each function in an environment of file, and 
#' identifies problem files, in order to help you refactor your code.
#' If the input is a path to an R file, then that file is sourced into 
#' a new environment and and the report is generated from that.
#' The number of lines of code that a function takes up is subjective 
#' in R; this function uses \code{length(deparse(fn))}. 
#' @examples
#' #Summarise function complexity in an environment
#' sig_report(pkg2env(stats))    
#' #Summarise function complexity in a file
#' \dontrun{
#' tmp <- tempfile(fileext = ".R")
#' writeLines(c(toString(sig(scan)), deparse(body(scan))), tmp)
#' sig_report(tmp)
#' }   
#' # Adjust the cutoff for reporting
#' sig_report(              
#'   baseenv(),  
#'   too_many_args  = 20,    
#'   too_many_lines = 100
#' )                
#' @export
sig_report <- function(x, ...)
{
  UseMethod("sig_report")
}
     
#' @rdname sig_report  
#' @method sig_report default
#' @export
sig_report.default <- function(x, ...)
{
  x <- as.environment(x)
  NextMethod("sig_report")
}

#' @rdname sig_report            
#' @method sig_report environment
#' @export
sig_report.environment <- function(x, too_many_args = 10, too_many_lines = 50, ...)
{
  all_vars <- mget(ls(envir = x, all.names = TRUE), envir = x)
  all_fns <- Filter(is.function, all_vars)
  n_args <- vapply(
    all_fns,
    function(fn) length(formals(fn)),
    integer(1)
  )
  n_lines <- vapply(
    all_fns,
    function(fn) length(deparse(fn)),
    integer(1)    
  )
  
  structure(
    list(
      n_vars             = length(all_vars),
      n_fns              = length(all_fns),
      n_args             = table(unname(n_args)),
      too_many_args      = too_many_args,
      fns_with_many_args = names(all_fns)[n_args > too_many_args],
      n_lines            = table(exponential_cut(n_lines)),
      too_many_lines     = too_many_lines,
      long_fns           = names(all_fns)[n_lines > too_many_lines]
    ),
    class = c("sigreport", "list")
  )   
}
                
#' @rdname sig_report           
#' @method sig_report character
#' @export
sig_report.character <- function(x, ...)
{  
  e <- source_to_new_env(x)
  sig_report(e, ...)
}

#' @rdname sig_report               
#' @method print sigreport
#' @export
print.sigreport <- function(x, ...)
{
  with(
    x, 
    cat(
      "The environment contains",
      n_vars,
      ngettext(n_vars, "variable", "variables"),
      "of which",
      n_fns,    
      ngettext(n_fns, "is a function.", "are functions."),
      "\nDistribution of the number of input arguments to the functions:"
    )
  )  
  print(x$n_args, ...)
  if(length(x$fns_with_many_args) > 0)
  {
    cat(
      "These functions have more than", 
      x$too_many_args, 
      "input args:\n"
    )
    print(noquote(x$fns_with_many_args), ...)
  } else
  {
    cat(
      "There are no functions with more than", 
      x$too_many_args, 
      "input args.\n"
    )
  }
  cat("Distribution of the number of lines of the functions:") 
  print(x$n_lines, ...)                                                 
  if(length(x$long_fns) > 0)
  {
    cat(
      "These functions have more than", 
      x$too_many_lines, 
      "lines:\n"
    )
    print(noquote(x$long_fns), ...)
  } else
  {
    cat(
      "There are no functions with more than", 
      x$too_many_lines, 
      "lines.\n"
    )
  }
}
