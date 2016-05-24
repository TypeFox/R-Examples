.BLEnv <- new.env()
.onLoad <- function(libname, pkgname)
{

  .BLEnv$settings <- list(gWidgetsToolkit = "tcltk", regFunc = "lm", 
            numSimulations = 50000, unitTestPath = system.file("RUnit", package = "BLCOP"))
}
