#' Functions to determine the operating system.
#'
#' \itemize{
#' \item{getOperatingSystem}{Simple wrapper for \code{.Platform$OS.type}, returns \code{character(1)}.}
#' \item{isUnix}{Predicate for OS string, returns \code{logical(1)}. Currently this would include Unix, Linux and Mac flavours.}
#' \item{isLinux}{Predicate for sysname string, returns \code{logical(1)}.}
#' \item{isDarwin}{Predicate for sysname string, returns \code{logical(1)}.}
#' \item{isWindows}{Predicate for OS string, returns \code{logical(1)}.}
#' }
#'
#' @return See above.
#' @export
getOperatingSystem = function() {
  .Platform$OS.type
}

#' @rdname getOperatingSystem
#' @export
isWindows = function() {
  .Platform$OS.type == "windows"
}

#' @rdname getOperatingSystem
#' @export
isUnix = function() {
  .Platform$OS.type == "unix"
}

#' @rdname getOperatingSystem
#' @export
isLinux = function() {
  isUnix() && grepl("linux", Sys.info()["sysname"], ignore.case = TRUE)
}

#' @rdname getOperatingSystem
#' @export
isDarwin = function() {
  isUnix() && grepl("darwin", Sys.info()["sysname"], ignore.case = TRUE)
}
