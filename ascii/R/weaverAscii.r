## Adapted from weaver package code by Seth Falcon

##' weaverAsciiFinish
##'
##' @param object object
##' @param error error
##' @keywords internal
weaverAsciiFinish <- function(object, error=FALSE) {
    resetStorage(cache_expr)
    RweaveAsciiFinish(object, error)
}

weaverAsciiRuncode <- makeRweaveAsciiCodeRunner(evalFunc=weaverEvalWithOpt)

##' weaverAsciidocSetup
##'
##' @param use.cache use.cache
##' @rdname RweaveAsciidoc
##' @keywords internal
weaverAsciiSetup <-
    function(file, syntax,
             output = NULL, quiet = FALSE, debug = FALSE,
             extension = "txt", backend = "asciidoc", openSchunk = "----",
             closeSchunk = "\n----\n", openSinput = "", closeSinput = "", openSoutput = "\n",
             closeSoutput = "", indent = "", openInclude = "include::", closeInclude = ".txt[]",
             openFig = "image::", closeFig = "[]", use.cache = TRUE, ...)
{
    if (!quiet)
      cat("Working dir:", getwd(), "\n")
    log_debug(paste("Working dir:", getwd()))
    res <- RweaveAsciiSetup(file = file, syntax = syntax, output = output, quiet = quiet, debug = debug, extension = extension, backend = backend, openSchunk = openSchunk, closeSchunk = closeSchunk, openSinput = openSinput, closeSinput = closeSinput, openSoutput = openSoutput, closeSoutput = closeSoutput, indent = indent, openInclude = openInclude, closeInclude = closeInclude, openFig = openFig, closeFig = closeFig, use.cache = use.cache, ...)
    res$options[["use.cache"]] <- use.cache
    res$options[["cache"]] <- FALSE
    ## to be on the safe side: see if defaults pass the check
    res$options <- RweaveAsciiOptions(res$options)
    res
}

##' weaverAsciidoc
##'
##' @rdname RweaveAsciidoc
##' @export
weaverAsciidoc <- function()
{
    require(weaver)
    list(setup = weaverAsciiSetup,
         runcode = weaverAsciiRuncode,
         writedoc = RweaveAsciiWritedoc,
         finish = weaverAsciiFinish,
         checkopts = RweaveAsciiOptions)
}

##' weaverT2tSetup
##'
##' @param use.cache use.cache
##' @rdname RweaveT2t
##' @keywords internal
weaverT2tSetup <- weaverAsciiSetup
formals(weaverT2tSetup) <- c(formals(RweaveT2tSetup)[1:(length(formals(RweaveT2tSetup))-1)], alist(use.cache=TRUE, ...=))

##' weaverT2t
##'
##' @rdname RweaveT2t
##' @export
weaverT2t <- function()
{
    require(weaver)
    list(setup = weaverT2tSetup,
         runcode = weaverAsciiRuncode,
         writedoc = RweaveAsciiWritedoc,
         finish = weaverAsciiFinish,
         checkopts = RweaveAsciiOptions)
}

##' weaverOrgSetup
##'
##' @param use.cache use.cache
##' @rdname RweaveOrg
##' @keywords internal
weaverOrgSetup <- weaverAsciiSetup
formals(weaverOrgSetup) <- c(formals(RweaveOrgSetup)[1:(length(formals(RweaveOrgSetup))-1)], alist(use.cache=TRUE, ...=))


##' weaverOrg
##'
##' @rdname RweaveOrg
##' @export
weaverOrg <- function()
{
    require(weaver)
    list(setup = weaverOrgSetup,
         runcode = weaverAsciiRuncode,
         writedoc = RweaveAsciiWritedoc,
         finish = weaverAsciiFinish,
         checkopts = RweaveAsciiOptions)
}

##' weaverPandocSetup
##'
##' @param use.cache use.cache
##' @rdname RweavePandoc
##' @keywords internal
weaverPandocSetup <- weaverAsciiSetup
formals(weaverPandocSetup) <- c(formals(RweavePandocSetup)[1:(length(formals(RweavePandocSetup))-1)], alist(use.cache=TRUE, ...=))

##' weaverPandoc
##'
##' @rdname RweavePandoc
##' @export
weaverPandoc <- function()
{
    require(weaver)
    list(setup = weaverPandocSetup,
         runcode = weaverAsciiRuncode,
         writedoc = RweaveAsciiWritedoc,
         finish = weaverAsciiFinish,
         checkopts = RweaveAsciiOptions)
}

##' weaverTextileSetup
##'
##' @param use.cache use.cache
##' @rdname RweaveTextile
##' @keywords internal
weaverTextileSetup <- weaverAsciiSetup
formals(weaverTextileSetup) <- c(formals(RweaveTextileSetup)[1:(length(formals(RweaveTextileSetup))-1)], alist(use.cache=TRUE, ...=))

##' weaverTextile
##'
##' @rdname RweaveTextile
##' @export
weaverTextile <- function()
{
    require(weaver)
    list(setup = weaverTextileSetup,
         runcode = weaverAsciiRuncode,
         writedoc = RweaveAsciiWritedoc,
         finish = weaverAsciiFinish,
         checkopts = RweaveAsciiOptions)
}

##' weaverReSTSetup
##'
##' @param use.cache use.cache
##' @rdname RweaveReST
##' @keywords internal
weaverReSTSetup <- weaverAsciiSetup
formals(weaverReSTSetup) <- c(formals(RweaveReSTSetup)[1:(length(formals(RweaveReSTSetup))-1)], alist(use.cache=TRUE, ...=))

##' weaverReST
##'
##' @rdname RweaveReST
##' @export
weaverReST <- function()
{
    require(weaver)
    list(setup = weaverReSTSetup,
         runcode = weaverAsciiRuncode,
         writedoc = RweaveAsciiWritedoc,
         finish = weaverAsciiFinish,
         checkopts = RweaveAsciiOptions)
}
