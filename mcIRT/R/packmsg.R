.onAttach <- function(libname, pkgname)
    {
      packageStartupMessage("Package: ",paste(pkgname)," is ready to estimate! \nThis is beta software! \nFollow this project on github: https://github.com/manuelreif/mcIRT.git")
    }