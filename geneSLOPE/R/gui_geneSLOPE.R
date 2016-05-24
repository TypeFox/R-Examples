#' GUI for GWAS with SLOPE
#'
#' A graphical user interface for performing Genome-wide
#' Association Study with SLOPE
#'
#' @details requires installing \pkg{\link[shiny]{shiny}} package
#'
#' @return null
#' @export
gui_geneSLOPE <- function() {
  appDir <- system.file("shiny-examples", "genSLOPE_gui", package = "geneSLOPE")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `geneSLOPE`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
