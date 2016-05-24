##' @export
`manifest` <-
function(x,...) UseMethod("manifest")

##' @export
`manifest.lvm` <-
function(x,...) {
  if (length(vars(x))>0)
    setdiff(vars(x),latent(x))
  else
    NULL
}

##' @export
`manifest.lvmfit` <-
function(x,...) {
  manifest(Model(x))
}

##' @export
manifest.list <- function(x,...) {
  manifestlist <- c()
  for (i in seq_along(x)) {
    manifestlist <- c(manifestlist, manifest(x[[i]]))
  }
  endolist <- unique(manifestlist)
  return(manifestlist)
}

##' @export
`manifest.multigroup` <-
function(x,...) {
  manifest(Model(x))
}
