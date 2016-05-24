".onAttach" <-
function (lib, pkg) 
{
  data(list=data(package=pkg)$results[,3],
       envir=as.environment(match(paste("package",pkg,sep=":"),
       search())))
 ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
 ver <- as.character(ver)
 packageStartupMessage("Rlab ", ver, " attached.\n")
} 
