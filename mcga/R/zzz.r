# Commented since R 2.15
#".First.lib" <-function (lib, pkg)
#{
#cat ("Loading package: Lib is ",lib, " and Pack is ",pack,"\n");
#Loading C library
# library.dynam ("mcga", pkg, lib);
#}


.onAttach <- function(lib, pkg){
  packageStartupMessage("Please use 'citation(\"mcga\")' for citing the R package mcga")
  invisible()
}

