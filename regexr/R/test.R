#' Test Regular Expression Validity
#' 
#' Test regular expression validity of a \code{regexr} object.
#' 
#' @param x A \code{regexr} object.
#' @param quiet logical.  Should \code{test} print warnings about the 
#' concatenated expression and individual sub-expressions?
#' @param \ldots Ignored.
#' @export
#' @return Returns a list of two logical vectors.  The first vector is a test of 
#' the concatenated expression.  The second vector is a logical test of the 
#' validity of each sub-expressions that makes up the concatenated 
#' expression.
#' @examples 
#' m <- construct(
#'     space = 
#'         "\\s+" 
#'             %:)%"I see",
#' 
#'     simp = 
#'         "(?<=(foo))",
#' 
#'     or = 
#'         "(;|:)\\s*"  
#'             %:)%"comment on what this does",
#' 
#'     "[a]s th[atey]"
#' )
#' 
#' 
#' test(m)
#' \dontrun{
#' subs(m)[5:7] <- c("(", "([A-Z]|(\\d{5})", ")")
#' test(m)
#' }
test <- function(x, quiet, ...){
    UseMethod("test")
}

