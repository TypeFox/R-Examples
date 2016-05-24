#' @name mgrepl
#' @keywords multiple grep sub gsub
#' @author Sven E. Templer
#' @title Multiple Pattern Matching and Replacement
#' @description 
#' \code{mgrepl} searches for any or all patterns and returns logical values.
#' Combination of the results is done via the logic functions \code{any},
#' \code{all} or \code{identity}, for example. 
#' Multicore feature is made available by \code{parallel::mclapply}.
#' @param patterns A character vector containing a regular expression
#' (\link{regex}) to be searched in \code{text}.
#' @param text Character vector where the search and replace is performed.
#' @param log.fun Logical function (\code{any} or \code{all}) 
#' to evaluate occurence of each pattern in \code{patterns} in each
#' value of \code{text}. Can also be custom. See examples.
#' @param use.which Logical, \code{TRUE} to return an integer like \link{which}
#' instead a logical vector.
#' @param cores Numeric value for how many cores to use for computation using 
#' \code{mclapply}.
#' @param \dots Further arguments passed to functions \code{grepl()}, \code{sub()} and \code{gsub()}.
#' @return
#' Logical vector of sam length as \code{text} where \code{TRUE} means either
#' any or all patterns in \code{patternlist} are matched in \code{text}
#' depending on \code{log.fun}.
#' Using \code{identity} as logical returns a matrix with results of patterns per row.
#' @seealso
#' \link{grep}, \link{mclapply}
#' @examples
#' #
#' 
#' # Compare different "log.fun" parameters:
#' mgrepl(c("a","b"), c("ab","ac","bc"), any)
#' mgrepl(c("a","b"), c("ab","ac","bc"), all)
#' mgrepl(c("a","b"), c("ab","ac","bc"), all, use.which = TRUE)
#' mgrepl(c("a","b"), c("ab","ac","bc"), identity)
#' mgrepl(letters[1:3], c("ax","xab","xbc"), function (x) sum(x)>1)
#' 
#' #

#' @export mgrepl
mgrepl <- function(patterns, text, log.fun = any, use.which = FALSE, cores = 1, ...) {

  ina <- is.na(text)
  patterns <- as.list(unlist(patterns))
  f <- match.fun(log.fun)
  i <- mclapply(patterns, function (y) grepl(y, text, ...), mc.cores=cores)
  i <- do.call(cbind, i)
  i <- apply(i, 1, f)
  if (use.which)
    return(which(i))
  i[ina] <- NA
  return(i)
  
}
