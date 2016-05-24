.onAttach=function(libname,pkgname){
   packageStartupMessage("Loaded gam ", as.character(packageDescription("gam")[["Version"]]),"\n")
}
