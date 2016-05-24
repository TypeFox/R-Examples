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
"LS2Wstat: Stationarity testing for LS2W fields\n\n",
"--- Written by Sarah Taylor and Matt Nunes ---\n",
"    --- Contributions from Idris Eckley ---\n",
"Current package version: ",ver," (",curdate,") \n\n",
"\n",
"**********************************************\n",
"\n",

"LS2Wstat", ver, "loaded\n")
)

}
