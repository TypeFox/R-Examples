#' @importFrom grDevices windowsFonts
.onAttach <- function(libname, pkgname) {

  if (!interactive()) return()
  Rcmdr <- options()$Rcmdr
  options(
    kmg2FontSize   = "14",
    kmg2FontFamily = 0,
    kmg2ColourSet  = 0,
    kmg2SaveGraph  = 0,
    kmg2Theme      = 0
  )
  
  if (.Platform$OS.type == "windows") {
    windowsFonts(
      Japan1 = "Japan1", Japan1HeiMin = "Japan1HeiMin",
      Japan1GothicBBB = "Japan1GothicBBB", Japan1Ryumin = "Japan1Ryumin",
      Korea1 = "Korea1", Korea1deb = "Korea1deb", CNS1 = "CNS1", GB1 = "GB1"
      )
  }
    
  plugins <- Rcmdr$plugins
  if (!pkgname %in% plugins) {
    Rcmdr$plugins <- c(plugins, pkgname)
    options(Rcmdr = Rcmdr)
    if("package:Rcmdr" %in% search()) {
      if(!getRcmdr("autoRestart")) {
        closeCommander(ask = FALSE, ask.save = TRUE)
        Commander()
      }
    }
    else {
      Commander()
    }
  }
  
}
