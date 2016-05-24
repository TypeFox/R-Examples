.onAttach <- function(libname,pkgname){
   packageStartupMessage("Loaded Iboot ", as.character(packageDescription("Iboot")[["Version"]]),"\n")
}

