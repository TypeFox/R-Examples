#
#  zzz.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2016 Bryon Aragam. All rights reserved.
#

# Export everything in this package so the other packages can use all of the methods
# #' @exportPattern ^[^\\.] <--- DON'T USE THIS!

#' @importFrom grDevices dev.off n2mfrow
#' @importFrom graphics par plot

.onAttach <- function(libname, pkgname){
    ### Only sparsebn needs a package startup message
    # packageStartupMessage("NOTE: This package is currently in a development state and may be unstable.\n Please report any bugs at https://github.com/itsrainingdata/sparsebnUtils/issues.")
}

### Set global options for this packages: Shamelessly stolen from /devtools/R/zzz.R
.onLoad <- function(libname, pkgname) {
  opt <- options()
#   op.devtools <- list(
#     devtools.path = "~/R-dev",
#     devtools.install.args = "",
#     devtools.name = "Your name goes here",
#     devtools.desc.author = 'person("First", "Last", email = "first.last@example.com", role = c("aut", "cre"))',
#     devtools.desc.license = "What license is it under?",
#     devtools.desc.suggests = NULL,
#     devtools.desc = list(),
#     devtools.revdep.libpath = file.path(tempdir(), "R-lib")
#   )

  opt.sparsebn <- list(
      sparsebn.graph = NULL,                              # graph package to use / NULL = no external package
      sparsebn.plotting = "igraph",                       # plotting mechanism to use
      sparsebn.zerothreshold = sqrt(.Machine$double.eps)  # threshold for numerical precision of zeroes
  )

  ### Only set the options which have not already been set (i.e. don't override defaults)
  toset <- !(names(opt.sparsebn) %in% names(opt)) # Which options have NOT been set already?
  if(any(toset)) options(opt.sparsebn[toset])     # Set the options which have not been set

  invisible()
}
