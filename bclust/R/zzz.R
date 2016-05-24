".onLoad" <-function(lib, pkg)
{
  library.dynam("bclust", package = pkg, lib.loc = lib)
  return(invisible(0)) 
}
plotnode<-getFromNamespace("plotNode","stats")
