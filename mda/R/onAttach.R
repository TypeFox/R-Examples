.onAttach=function(libname,pkgname){
   packageStartupMessage("Loaded mda ", as.character(packageDescription("mda")[["Version"]]),"\n")
}
