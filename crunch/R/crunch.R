##' Crunch.io: instant, visual, collaborative data analysis
##'
##' \href{http://crunch.io/}{Crunch.io} provides a cloud-based data store and
##' analytic engine. It has a \href{https://beta.crunch.io/}{web client} for
##' interactive data exploration and visualization. The crunch package for R
##' allows analysts to interact with and manipulate Crunch datasets from within
##' R. Importantly, this allows technical researchers to collaborate naturally
##' with team members, managers, and clients who prefer a point-and-click
##' interface: because all connect to the same dataset in the cloud, there is no
##' need to email files back and forth continually to share results.
##'
##' @seealso To learn more about using the package, see
##' \code{vignette("getting-started", package="crunch")}. To sign up for a
##' Crunch.io account, visit \url{https://beta.crunch.io/}.
##' @docType package
##' @name crunch
NULL

##' @importFrom httr set_config
.onAttach <- function (lib, pkgname="crunch") {
    setIfNotAlready(
        crunch.api="https://beta.crunch.io/api/",
        crunch.max.categories=256,
        crunch.timeout=60,
        httpcache.on=TRUE,
        crunch.namekey.dataset="alias",
        crunch.namekey.array="name"
    )
    set_config(crunchConfig())
    notifyIfNewVersion()
    invisible()
}

setIfNotAlready <- function (...) {
    newopts <- list(...)
    oldopts <- options()
    oldopts <- oldopts[intersect(names(newopts), names(oldopts))]
    newopts <- updateList(newopts, oldopts)
    do.call(options, newopts)
    invisible(oldopts)
}
