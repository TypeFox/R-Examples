#' @name match_arg
#' @export
#' 
#' @title Argument Verification Using Partial Matching
#' @description \code{match_arg} matches \code{arg} against a table of candidate values as 
#'   specified by \code{choices}, where \code{NULL} means to take the first one.  This is a modified
#'   version of the \code{base} package function \code{match.arg} that differs only in that it 
#'   can interact with \code{ArgumentCheck} environments.
#'   
#' @param arg a character vector (of length one unless \code{several.ok} is \code{TRUE}) or 
#'   \code{NULL}.
#' @param choices a character vector of candidate values
#' @param several.ok logical specifying if \code{arg} should be allowed to have more than one element.
#' @param argcheck An \code{ArgumentCheck} environment as returned by \code{\link{newArgCheck}}
#' 
#' @details In the one-argument form \code{match_arg(arg)}, the \code{choices} are obtained from 
#'   a default setting for the formal argument \code{arg} of the function from which 
#'   \code{match.arg} was called. (Since default argument matching will set \code{arg} to 
#'   \code{choices}, this is allowed as an exception to the 'length one unless \code{several.ok} 
#'   is \code{TRUE}' rule, and returns the first element.)
#'   
#'   Matching is done using \code{\link{pmatch}}, so \code{arg} may be abbreviated.
#'   
#' @return The unabbreviated version of the exact or unique partial match if there is one; 
#'   otherwise, an error is signalled if \code{several.ok} is false, as per default. 
#'   When \code{several.ok} is true and more than one element of \code{arg} has a match, 
#'   all unabbreviated versions of matches are returned.
#'   
#' @author R Core.  This function is a near-verbatim copy of the \code{\link{match.arg}} in 
#'   \code{base} R.
#'   
#' @seealso 
#' \code{\link{pmatch}}, \code{\link{match.fun}}, \code{\link{match.call}}
#' 
#' @examples 
#' require(stats)
#' # Extends the example for 'switch'
#' center <- function(x, type = c("mean", "median", "trimmed")) {
#'   type <- match.arg(type)
#'     switch(type,
#'            mean = mean(x),
#'            median = median(x),
#'            trimmed = mean(x, trim = .1))
#' }
#' x <- rcauchy(10)
#' center(x, "t")       # Works
#' center(x, "med")     # Works
#' try(center(x, "m"))  # Error
#' stopifnot(identical(center(x),       center(x, "mean")),
#'           identical(center(x, NULL), center(x, "mean")) )
#'
#' ## Allowing more than one match:
#' match.arg(c("gauss", "rect", "ep"),
#'           c("gaussian", "epanechnikov", "rectangular", "triangular"),
#'           several.ok = TRUE)


match_arg <- function (arg, choices, several.ok = FALSE, argcheck) 
{
  if (missing(choices)) {
    formal.args <- formals(sys.function(sys.parent()))
    choices <- eval(formal.args[[deparse(substitute(arg))]])
  }
  if (is.null(arg)) 
    return(choices[1L])
  else if (!is.character(arg)) 
    addError(paste0("'", deparse(substitute(arg)), "' must be NULL or a character vector"),
             argcheck)
  if (!several.ok) {
    if (identical(arg, choices)) 
      return(arg[1L])
    if (length(arg) > 1L) 
      addError(paste0("'", deparse(substitute(arg)), "' must be of length 1"),
               argcheck)
  }
  else if (length(arg) == 0L) 
    addError(paste0("'", substitute(arg), "' must be of length >= 1"),
             argcheck)
  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L)) 
    addError(gettextf(paste0("'", deparse(substitute(arg)), "' should be one of %s"), 
                      paste(dQuote(choices),
                            collapse = ", "),
                      domain = NA),
             argcheck)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1) 
    addError("there is more than one match in 'match_arg'",
             argcheck)
  choices[i]
}