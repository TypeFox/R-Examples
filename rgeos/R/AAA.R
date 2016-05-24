.RGEOS_HANDLE <- new.env(FALSE, parent=globalenv())

set_RGEOS_HANDLE <- function(handle) {
    assign("GEOSptr", handle, envir=.RGEOS_HANDLE)
}

set_RGEOS_DENSE <- function(value) {
    stopifnot(is.logical(value))
    stopifnot(length(value) == 1)
    assign("returnDense", value, envir=.RGEOS_HANDLE)
}

get_RGEOS_DENSE <- function() {
    get("returnDense", envir=.RGEOS_HANDLE)
}

set_RGEOS_dropSlivers <- function(value) {
    stopifnot(is.logical(value))
    stopifnot(length(value) == 1)
    assign("dropSlivers", value, envir=.RGEOS_HANDLE)
}

get_RGEOS_dropSlivers <- function() {
    get("dropSlivers", envir=.RGEOS_HANDLE)
}

set_RGEOS_warnSlivers <- function(value) {
    stopifnot(is.logical(value))
    stopifnot(length(value) == 1)
    assign("warnSlivers", value, envir=.RGEOS_HANDLE)
}

get_RGEOS_warnSlivers <- function() {
    get("warnSlivers", envir=.RGEOS_HANDLE)
}

set_RGEOS_polyThreshold <- function(value) {
    stopifnot(is.numeric(value))
    stopifnot(length(value) == 1)
    stopifnot(value >= 0.0)
    assign("polyThreshold", value, envir=.RGEOS_HANDLE)
}

get_RGEOS_polyThreshold <- function() {
    get("polyThreshold", envir=.RGEOS_HANDLE)
}

set_RGEOS_STR <- function(value) {
    stopifnot(is.logical(value))
    stopifnot(length(value) == 1)
    assign("STRsubset", value, envir=.RGEOS_HANDLE)
}

get_RGEOS_STR <- function() {
    get("STRsubset", envir=.RGEOS_HANDLE)
}

init_RGEOS <- function() {
    .Call('rgeos_Init', PACKAGE="rgeos")
}

finish_RGEOS <- function() {
    .Call('rgeos_finish', .RGEOS_HANDLE, PACKAGE="rgeos")
}

version_GEOS <- function() {
    .Call("rgeos_GEOSversion", PACKAGE="rgeos")
}

version_GEOS0 <- function() {
    substring(version_GEOS(), 1, 5)
}

version_sp_linkingTo <- function() {
    .Call("rgeos_sp_linkingTo_version")
}

.onLoad <- function(lib, pkg) {
#  require(methods, quietly = TRUE, warn.conflicts = FALSE)
#  require("sp")
#  require("stringr")
#  library.dynam('rgeos', pkg, lib)

  set_RGEOS_HANDLE(init_RGEOS())
  assign("scale", 100000000, envir=.RGEOS_HANDLE)
  assign("do_poly_check", TRUE, envir=.RGEOS_HANDLE)
#  assign("both_poly", FALSE, envir=.RGEOS_HANDLE)
#  assign("drop_not_poly", FALSE, envir=.RGEOS_HANDLE)
  assign("polyThreshold", 0.0, envir=.RGEOS_HANDLE)
  assign("dropSlivers", FALSE, envir=.RGEOS_HANDLE)
  assign("warnSlivers", TRUE, envir=.RGEOS_HANDLE)
  assign("returnDense", TRUE, envir=.RGEOS_HANDLE)
  assign("STRsubset", FALSE, envir=.RGEOS_HANDLE)
}

.onAttach <- function(lib, pkg) {
  fn <- system.file("SVN_VERSION", package="rgeos")
  if (file.exists(fn)) {
    svn_version <- scan(fn, what=character(1), sep="\n", quiet=TRUE)
  } else {
    svn_version <- "(unknown)"
  }
  Smess <- paste("rgeos version: ", utils::packageDescription("rgeos")$Version,
    ", (SVN revision ", svn_version, ")\n", sep="")
  Smess <- paste(Smess, "GEOS runtime version:",
    version_GEOS(), "\n")
  splVersion <- version_sp_linkingTo()
  Smess <- paste(Smess, "Linking to sp version:", splVersion, "\n")
  spVcheck <- NULL
  if("sp" %in% .packages()) spVcheck <- utils::packageVersion("sp") == splVersion
  if (!is.null(spVcheck) && !spVcheck) paste(Smess, 
    "sp version used to install rgeos and loaded sp version differ\n")

  Smess <- paste(Smess, "Polygon checking:", get_do_poly_check(), "\n")
  packageStartupMessage(Smess, appendLF = TRUE)
}

.onUnload <- function(libpath) {
  invisible(finish_RGEOS())
}
