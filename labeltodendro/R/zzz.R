".onLoad" <-function(lib, pkg)
{
  library.dynam("labeltodendro", package = pkg, lib.loc = lib)
  return(invisible(0)) 
}
