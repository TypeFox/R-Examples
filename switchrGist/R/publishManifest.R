##' @import httpuv
##' @import gistr
##' @import switchr
##' @import methods
NULL

##' publishManifest
##'
##' publish a manifest to Github in the form of a gist
##'
##' @param manifest The manifest to be published
##' @param dest a GistDest object
##' @param desc description to apply to the Gist
##' @param fname The name of the file to create within the gist.
##' @param ... unused
##' @return the url to access the raw file within the gist.
##' @docType methods
##' @aliases publishManifest,PkgManifest,GistDest
##' @export
##' @examples
##' man = libManifest()
##'\dontrun{gisturl = publishManifest(man, Gist())}
##' @rdname publishManifest
##' @importFrom RJSONIO toJSON
setMethod("publishManifest", c(manifest = "PkgManifest",
                               dest = "GistDest"),
          function(manifest, dest, desc = "An R package manifest", fname = "manifest.rman",...) {
              .publishManifest2(man = manifest,
                               desc = desc)
          })
##' @export
##' @rdname publishManifest
##' @aliases publishManifest,SessionManifest,GistDest

setMethod("publishManifest", c(manifest = "SessionManifest",
                               dest = "GistDest"),
          function(manifest, dest, desc="An R seeding manifest", fname = "manifest.rman", ...) {
              .publishManifest2(man = manifest,
                               desc = desc)
          })

.publishManifest2 = function(man, desc = "An R package manifest", fname = "manifest.rman", ...) {
    fil = tempfile(pattern = fname)
    fil = publishManifest(man, fil)
    g = gist_create(files = fil, description = desc)
    fname = names(g$files)
    g$files[[fname]]$raw_url
}

    



