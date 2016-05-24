##' Set or return options of phylobase
##' 
##' Provides a mean to control the validity of \code{phylobase}
##' objects such as singletons, reticulated trees, polytomies, etc.
##' 
##' The parameter values set via a call to this function will remain
##' in effect for the rest of the session, affecting the subsequent
##' behavior of phylobase.
##' 
##' @param \dots a list may be given as the only argument, or any
##' number of arguments may be in the \code{name=value} form, or no
##' argument at all may be given.  See the Value and Details sections
##' for explanation.
##' @return A list with the updated values of the parameters. If
##' arguments are provided, the returned list is invisible.
##' @author Francois Michonneau (adapted from the package \code{sm})
##' @keywords phylobase validator
##' @examples
##' \dontrun{
##' phylobase.options(poly="fail")
##' # subsequent trees with polytomies will fail the validity check
##' }
##' 
##' @export
phylobase.options <- function (...) {
    if (nargs() == 0) return(.phylobase.Options)
    current <- .phylobase.Options
    temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list = temp <- arg,
               character = return(.phylobase.Options[arg]),
               stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")

    if (!all(names(temp) %in% names(current)))
        stop("Option name invalid: ", sQuote(names(temp)))
    changed <- current[n]
    current[n] <- temp
    current <- lapply(current, function(foo) {
        foo <- match.arg(foo, c("warn", "fail", "ok"))
    })
    if (!identical(current$retic, "fail")) {
        stop("Currently reticulated trees are not handled by phylobase.")
    }
    ## options are always global
    env <- asNamespace("phylobase")
    assign(".phylobase.Options", current, envir = env)
    invisible(current)
}
