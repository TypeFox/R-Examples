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
"nlt: a package to perform nondecimated wavelet lifting\n\n",
"--- Written by Marina Knight and Matt Nunes ---\n",
"  Current package version: ",ver," (",curdate,") \n\n",
"            -+ packaged by MAN +-           \n",
"**********************************************\n",
"\n",

"nlt", ver, "loaded\n")
)

}
