## Show errata

##' Show errata
##'
##' @return opens browse to errata page
##' @export
errata <- function() {
  f <- system.file("errata", "errata.html", package="UsingR")
  browseURL(f)
}
