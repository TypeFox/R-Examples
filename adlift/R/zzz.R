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
"adlift: a package to perform wavelet lifting schemes\n\n",
"--- Written by Matt Nunes and Marina Knight ---\n",
"  Current package version: ",ver," (",curdate,") \n\n",
"            -+ packaged by MAN +-           \n",
"**********************************************\n",
"\n",

"adlift", ver, "loaded\n")
)

}
