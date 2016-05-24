.onAttach <-
function (libname, pkgname) 
{
  msg <- paste("Loading 'metasens' package (version ",
               utils::packageDescription("metasens")$Version,
               ").", sep="")
  packageStartupMessage(msg)
}
