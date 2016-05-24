.onUnload <- function(libpath)
{
  if (is.loaded("treatSens_fitSensitivityAnalysis", PACKAGE = "treatSens")) {
    library.dynam.unload("treatSens", libpath)
  }
}
