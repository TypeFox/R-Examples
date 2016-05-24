.onAttach <- function(libname, pkgname){
    packageStartupMessage("tuneR >= 1.0 has changed its Wave class definition.\nUse updateWave(object) to convert Wave objects saved with previous versions of tuneR.")
}
