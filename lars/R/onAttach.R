.onAttach=function(libname,pkgname){
   packageStartupMessage("Loaded lars ", as.character(packageDescription("lars")[["Version"]]),"\n")
}
