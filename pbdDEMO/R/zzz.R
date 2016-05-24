.onLoad <- function(libname, pkgname)
{
  pbdMPI::pbd_opt(gbd.major = 1L)
  pbdMPI::pbd_opt(divide.method = c("block.cyclic", "block0"))
  
  invisible()
}

