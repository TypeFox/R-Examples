# add libraries to ressource path
.onLoad <- function(libname, pkgname) {
  shiny::addResourcePath("RLumShiny", system.file("www", package = "RLumShiny"))
}

# Dependencies in the shiny apps are currently not registered by R CMD check --as-cran
.satisfyCheck <- function() {
  x <- TRUE
  if (x) return(x)
  Luminescence::sTeve()
  RCurl::reset()
  googleVis::renderGvis()
  digest::digest(".")
}