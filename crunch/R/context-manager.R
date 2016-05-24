##' Context managers
##'
##' @param enter function to run before doing things
##' @param exit function to run after doing things
##' @param error optional function to run if an error is thrown
##' @param as character optional way to specify a default name for assinging
##' the return of the enter function.
##' @return an S3 class "contextManager" object
##' @seealso \link{with-context-manager}
##' @aliases contextManager
##' @export
ContextManager <- function (enter=function (){}, exit=function (){},
                            error=NULL, as=NULL) {
    structure(list(enter=enter, exit=exit, error=error, as=as),
        class="contextManager")
}

##' Context manager's "with" method
##'
##' @param data \code{\link{contextManager}}
##' @param expr code to evaluate within that context
##' @param ... additional arguments. One additional supported argument is "as",
##' which lets you assign the return of your "enter" function to an object you
##' can access.
##' @return Nothing.
##' @name with-context-manager
##' @seealso \link{ContextManager}
##' @export
with.contextManager <- function (data, expr, ...) {
    env <- parent.frame()
    on.exit(data$exit())
    setup <- data$enter()
    dots <- list(...)
    as.name <- dots$as %||% data$as
    if (!is.null(as.name)) {
        assign(as.name, setup, envir=env)
        ## rm this after running? or add the rm step to the exit
    }
    if (is.function(data$error)) {
        tryCatch(eval(substitute(expr), envir=parent.frame()), error=data$error)
    } else {
        eval(substitute(expr), envir=parent.frame())
    }
}

##' Set some global options temporarily
##'
##' @param ... named options to set
##' @return an S3 class "contextManager" object
##' @seealso \link{with-context-manager} \link{ContextManager}
##' @export
temp.options <- function (...) {
    new <- list(...)
    old <- sapply(names(new), getOption, simplify=FALSE)
    return(ContextManager(
        function () do.call(options, new),
        function () do.call(options, old)
    ))
}

##' @rdname temp.options
##' @export
temp.option <- temp.options

##' Give consent to do things that require permission
##' @return an S3 class "contextManager" object
##' @seealso \link{with-context-manager} \link{ContextManager}
##' @export
consent <- function () {
    temp.options(crunch.require.confirmation=FALSE)
}

requireConsent <- function () {
    getOption("crunch.require.confirmation") %||% interactive()
}
