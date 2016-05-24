.onAttach <-
function (libname, pkgname) 
{
  msg <- paste("Loading 'netmeta' package (version ",
               utils::packageDescription("netmeta")$Version,
               ").", sep="")
  packageStartupMessage(msg)
}
