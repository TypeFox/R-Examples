#' @import stats
#' @import checkmate

.onLoad = function(libname, pkgname) {
  options(BBmisc.ProgressBar.stream = getOption("BBmisc.ProgressBar.stream", "stderr"))
  options(BBmisc.ProgressBar.style = getOption("BBmisc.ProgressBar.style", "text"))
  options(BBmisc.ProgressBar.width = getOption("BBmisc.ProgressBar.width", getOption("width")))
}
