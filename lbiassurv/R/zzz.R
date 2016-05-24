".onLoad" <-function(lib, pkg)
{
  library.dynam("lbiassurv", package = pkg, lib.loc = lib)
  return(invisible(0)) 
}
