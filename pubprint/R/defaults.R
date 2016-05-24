#############################################################################
# defaults.R
#############################################################################

#' @include out.R
#' @include style.R
NULL

new_defaults <- function(value = list()) 
{
    defaults <- value

    get <- function(name) 
    {
        if (missing(name)) 
            defaults 
        else 
        {
            if (1 == length(name)) 
                defaults[[name]] 
            else
                setNames(defaults[name], name)
        }
    }

    set <- function(...)
    {
        dots = list(...)

        if (0 == length(dots)) return()
        if (is.null(names(dots)) && length(dots) == 1 && is.list(dots[[1]]))
            if (length(dots <- dots[[1]]) == 0) return()

        defaults[names(dots)] <<- dots
        invisible(NULL)
    }

    list(get = get, set = set)
}

#' Output format options for the pubprint package
#'
#' A list which functions are used to print in the correct output format
#' (LaTeX, HTML, Markdown or plain text).
#'
#' Using \code{pp_opts_out$get()} shows all currently used output format
#' functions, \code{pp_opts_out$set()} allows to change them.
#'
#' @seealso See \code{\link{pp_init_out}} for initialising this variable in the
#' correct way and \code{\link{pp_init_style}} for publication style.
#' 
#' @examples
#' pp_opts_out$set(pp_init_out())
#' pp_opts_out$set(pp_init_out("html"))
#' 
#' @export
pp_opts_out <- new_defaults(pp_init_out())

#' Publication style options for the pubprint package
#'
#' A list which functions are used to print in the correct publication style
#' (like APA).
#'
#' Using \code{pp_opts_style$get()} to show all currently used publication
#' style functions. \code{pp_opts_style$set()} allows to change them.
#'
#' @seealso See \code{\link{pp_init_style}} for initialising this variable in the
#' correct way and \code{\link{pp_init_out}} for the output format.
#' 
#' @examples
#' pp_opts_style$set(pp_init_style())
#' pp_opts_style$set(pp_init_style("apa"))
#' 
#' @export
pp_opts_style <- new_defaults(pp_init_style())

#' General options for the pubprint package
#'
#' Options including how many decimal places are used, whether to remove items
#' when pulling them, etc.
#'
#' Set global options with \code{pp_opts$set()} and get your options with
#' \code{pp_opts$get()}.
#' \describe{
#'   \item{\code{nsmall}:}{controls the number of digits to print when printing
#'     numeric values.}
#'   \item{\code{leading0}:}{controls whether a leading zero is printed. In
#'     some cases this may be inappropriate, like for correlation coefficients
#'     (r = .73).}
#'   \item{\code{replace0}:}{controls whether numbers that are absolute smaller than
#'     \code{1/(10^nsmall)} are replaced with this term. 
#'      For example \eqn{p = .0001} with \eqn{p < .001}. Pay attention that
#'      zero is replaced as well.}
#'   \item{\code{drop0trailing}:}{logical, indicating if trailing zeros, i.e.,
#'      \code{"0"} _after_ the decimal mark, should be removed; also drops
#'      \code{"e+00"} in exponential formats.}
#'   \item{\code{delimiter}:}{delimiter between items.}
#'   \item{\code{removeItems}:}{controls whether items are removed when pulling
#'     them. Either a logical, \code{"memory"} or \code{"pipe"}. See
#'     \code{\link{pull.pubprint} for more details.}}
#'   \item{\code{mmode}:}{controls whether output is set in math mode.}
#'   \item{\code{brackets}:}{controls which brackets are used.}
#'   \item{\code{separator}:}{controls whether a separator between content and
#'     pubprint output is printed.}
#' }
#'
#' @examples
#' pp_opts$set(nsmall = 3)
#' pp_opts$set(nsmall = 3, removeItems = FALSE)
#'
#' @export
pp_opts <- new_defaults(list(nsmall = 2,
                             leading0 = TRUE,
                             replace0 = TRUE,
                             drop0trailing = FALSE,
                             delimiter = ", ",
                             removeItems = "pipe",
                             mmode = TRUE,
                             brackets = c("(", ")", "[", "]"),
                             separator = "brackets"))
