#' Plot an FPCA object. 
#'
#' Currently implemented as just plotting the diagnostics plots, the same as CreateDiagnosticsPlot.
#' @param x An FPCA class object returned by FPCA()
#' @param ... passed into CreateDiagnosticsPlot
#' @export
#' @seealso CreateDiagnosticsPlot
plot.FPCA <- function(x, ...) {
  CreateDiagnosticsPlot(x, ...)
}