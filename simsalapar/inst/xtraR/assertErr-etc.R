if( exists("assertCondition", asNamespace("tools")) ) { ## R > 3.0.1

assertError <- function(expr, verbose=getOption("verbose"))
    tools::assertCondition(expr, "error", verbose=verbose)
assertWarning <- function(expr, verbose=getOption("verbose"))
    tools::assertCondition(expr, "warning", verbose=verbose)
assertWarningAtLeast <- function(expr, verbose=getOption("verbose"))
    tools::assertCondition(expr, "error", "warning", verbose=verbose)

} else { ## in R <= 3.0.1, use our old versions

##' @title Ensure evaluating 'expr' signals an error
##' @param expr
##' @return the caught error, invisibly
##' @author Martin Maechler
assertError <- function(expr, verbose=getOption("verbose")) {
    d.expr <- deparse(substitute(expr))
    t.res <- tryCatch(expr, error = function(e) e)
    if(!inherits(t.res, "error"))
	stop(d.expr, "\n\t did not give an error", call. = FALSE)
    cat("Asserted Error:", conditionMessage(t.res),"\n")
    invisible(t.res)
}

##' @title Ensure evaluating 'expr' signals a warning
##' @param expr
##' @return the caught warning, invisibly
##' @author Martin Maechler
assertWarning <- function(expr, verbose=getOption("verbose")) {
    d.expr <- deparse(substitute(expr))
    t.res <- tryCatch(expr, warning = function(w)w)
    if(!inherits(t.res, "warning"))
	stop(d.expr, "\n\t did not give a warning", call. = FALSE)
    if(verbose) cat("Asserted Warning:", conditionMessage(t.res),"\n")
    invisible(t.res)
}

}## if .. else { not yet in tools }
