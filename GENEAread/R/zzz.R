.onAttach <-function(lib,pkg)
{
ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
     ver <- as.character(ver)	

options(digits=12)
text = paste("GENEAread", ver, "loaded\n")
packageStartupMessage(text)

}
