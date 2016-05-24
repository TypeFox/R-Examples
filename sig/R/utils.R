#' Get environment of a package.
#' 
#' Utility function to get the environment of a package on the search 
#' path.
#' 
#' @param pkg A package.
#' @return the environment corresponding to \code{pkg}.
#' @seealso \code{\link[base]{list2env}}
#' @examples
#' pkg2env(graphics)
#' @export
pkg2env <- function(pkg) 
{
  pkg_name <- deparse(substitute(pkg))
  if(!pkg_name %in% .packages())
  {
    if(pkg_name %in% .packages(TRUE))
    {
      message("Loading package ", sQuote(pkg_name), ".")
      library(pkg_name, character.only = TRUE)
    } else
    {
      stop("The package ", sQuote(pkg_name), " is not available.")
    }
  }
  as.environment(paste0("package:", pkg_name))
}

#' Source a file into a new environment.
#' 
#' Silently sources a file into a new environment, 
#' returning that environment.
#' @param file a file to source.
#' @param encoding character encoding of that file.
#' @return An environment containing the sourced variables.
source_to_new_env <- function(file, encoding = getOption("encoding"))
{  
  e <- new.env()
  source(file, e, verbose = FALSE, encoding = encoding)
  e
}

#' Wrap in backquotes
#' 
#' Wraps strings in backquotes.
#' @param x A character vector.
#' @return A character vector.
#' @note Existing backquote characters are escaped with a backslash.
#' @seealso \code{\link[base]{sQuote}}
#' @examples
#' \dontrun{
#' backquote(c("foo bar", "a`b`c"))
#' }
backquote <- function(x)
{
  x <- gsub("`", "\\\\`", x)
  paste0("`", x, "`")
}

#' Fix names for sigs
#' 
#' Make anonymous functions and special functions safe.
#' @param fn_name A character vector.
#' @return A character vector.
#' @note Strings beginning with ``function'' are given the value
#' \code{"..anonymous.."}.
#' 
#' Special function names are wrapped in backquotes.
#' @examples
#' \dontrun{
#' fix_fn_names(c("%foo%", "?", "foo bar", "repeat", "function"))
#' }
fix_fn_names <- function(fn_name)
{
  fn_name[grepl("^function\\(", fn_name)] <- "..anonymous.."
  is_special <- make.names(fn_name) != fn_name
  fn_name[is_special] <- backquote(fn_name[is_special])
  fn_name
}

#' Cut with exponential breaks
#' 
#' Wrapper to \code{cut} for positive integers.
#' @param x A vector of positive integers.
#' @return A factor.
#' @note The breaks are 1, 2, 3 to 4, 5 to 8, etc. 
#' No input checking is done; use at your peril.
#' @seealso \code{\link[base]{cut}}
#' @examples
#' \dontrun{
#' exponential_cut(c(1:10, 500))
#' }
exponential_cut <- function(x)
{
  cut_points <- c(0, 2 ^ seq.int(0, ceiling(log2(max(x)))))
  n_cut_points <- length(cut_points)
  lo <- cut_points[-n_cut_points] + 1
  hi <- cut_points[-1]
  labels <- ifelse(
    lo == hi,
    lo,
    paste0("[", lo, ",", hi, "]")
  )
  cut(x, cut_points, labels = labels)
}

#' Workhorse of the print methods
#' 
#' Wraps toString methods with cat.
#' @param x Object to print
#' @param ... Passed to \code{toString}.
#' @return The input is invisibly returned, but the function is mostly invoked for the side effect of printing the object.
#' @note Not intended for general consumption.  This function is only 
#' exported because of package build requirements.
#' @export
print_engine <- function(x, ...)
{
  cat(toString(x, ...), sep = "\n")
  invisible(x)
}
