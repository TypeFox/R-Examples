#' Extended Argument Matching
#' 
#' \code{match.arg.ext} matches \code{arg} against a set of candidate values as specified by \code{choices}; it extends \code{\link{match.arg}} by allowing \code{arg} to be a numeric identifier of the \code{choices}.
#' @name match.arg.ext
#' @param arg a character string or numeric value
#' @param choices a character vector of candidate values
#' @param base numeric value, specifying the numeric index assigned to the first element of \code{choices}
#' @param several.ok logical specifying if \code{arg} should be allowed to have more than one element
#' @param numeric logical specifying if the function should return the numerical index (counting from \code{base}) of the matched \code{arg}ument, or, by default, its name 
#' @param ignore.case logical specifying if the matching should be case sensitive
#' @details When \code{choices} are missing, they are obtained from a default setting for the formal argument \code{arg} of the function from which \code{match.arg.ext} was called.
#'
#' Matching is done using \code{\link{pmatch}} (indirectly through a call to \code{\link{match.arg}}, so \code{arg} may be abbreviated.
#'
#' If \code{arg} is numeric, it may take values between \code{base} and \code{length(choices)+base-1}.  \code{base=1} will give standard 1-based R indices, \code{base=0} will give indices counted from zero as used to identify SAGA modules in library RSAGA.
#' @return If \code{numeric} is false and \code{arg} is a character string, the function returns the unabbreviated version of the unique partial match of \code{arg} if there is one; otherwise, an error is signalled if \code{several.ok} is false, as per default. When \code{several.ok} is true and there is more than one match, all unabbreviated versions of matches are returned.
#'
#' If \code{numeric} is false but \code{arg} is numeric, \code{match.arg.ext} returns name of the match corresponding to this index, counting from \code{base}; i.e. \code{arg=base} corresponds to \code{choices[1]}.
#'
#' If \code{numeric} is true, the function returns the numeric index(es) of the partial match of \code{arg}, counted from \code{base} to \code{length(choices)+base-1}. If \code{arg} is already numeric, the function only checks whether it falls into the valid range from \code{arg} to \code{length(choices)+base-1} and returns \code{arg}.
#' @author Alexander Brenning
#' @seealso \code{\link{match.arg}}, \code{\link{pmatch}}
#' @examples
#' # Based on example from 'match.arg':
#' require(stats)
#' center <- function(x, type = c("mean", "median", "trimmed")) {
#'   type <- match.arg.ext(type,base=0)
#'   switch(type,
#'          mean = mean(x),
#'          median = median(x),
#'          trimmed = mean(x, trim = .1))
#' }
#' x <- rcauchy(10)
#' center(x, "t")       # Works
#' center(x, 2)         # Same, for base=0
#' center(x, "med")     # Works
#' center(x, 1)         # Same, for base=0
#' try(center(x, "m"))  # Error
#' @keywords utilities
#' @export
match.arg.ext = function(arg, choices, base = 1, several.ok = FALSE, 
    numeric = FALSE, ignore.case = FALSE) 
{
    if (missing(choices)) {
        formal.args <- formals(sys.function(sys.parent()))
        choices <- eval(formal.args[[deparse(substitute(arg))]])
    }
    if (is.character(arg)) {
        if (ignore.case) {
            choices = tolower(choices)
            arg = tolower(arg)
        }
        res = match.arg(arg=arg,choices=choices,several.ok=several.ok)
        if (numeric)  res = which(choices %in% res) + base - 1
    } else if (is.numeric(arg)) {
        if ( (arg<base) | (arg>(length(choices)+base-1)) )
            stop("'arg' should be between ",base," and ",length(choices)+base-1)
        if (numeric) {
            res = arg
        } else {
            res = choices[arg - base + 1]
        }
    } else stop("'arg' should be numeric or character")
    return(res)
}
