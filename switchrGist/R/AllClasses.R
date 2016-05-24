
##'GistDest
##'
##' A simple object indicating that the publishing destination
##' is a GitHub gist.
##' @export
setClass("GistDest", representation(unused="logical"))

##' Gist
##'
##' Create a GistDest object. Pass this to publishManifest as
##' the destination to publish a package or seeding manifest to a
##' GitHub gist.
##' @export
Gist = function() new("GistDest")

setMethod("show", "GistDest",
          function(object) cat("\nA GistDest object.\n"))
