
.onAttach <- function(libname, pkgname)
{
  packageStartupMessage(paste(pkgname)," package calling ...\nFollow this project on github: https://github.com/manuelreif/PP.git")
}
