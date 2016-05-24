.onAttach <-function(lib,pkg)
{
ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
     ver <- as.character(ver)	

curdate <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Date")
    curdate <- as.character(curdate)

# Welcome message (MAN):

packageStartupMessage(paste(
"\n",
"**********************************************\n",
"binhf: Haar-Fisz functions for binomial data\n\n",
"--- Written by Matt Nunes ---\n",
"  Current package version: ",ver," (",curdate,") \n\n",
"\n",
"**********************************************\n",
"\n",

"binhf", ver, "loaded\n")
)

}
