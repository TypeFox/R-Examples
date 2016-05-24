
.onLoad <- function(lib, pkg){

  ##library.dynam("gRain", package = pkg, lib.loc = lib)
  return(invisible(0))
}

## .onAttach <- function(lib, pkg){

##   if('Rgraphviz' %in% rownames(installed.packages())){
##     library(Rgraphviz)
##   } else {
##     cat("Note: To display models the Rgraphviz package (from Bioconductor) must be installed.\n") 
##   }
  
## }
