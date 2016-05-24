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
"abctools: A package with tools for ABC inference\n\n",
"--- Written by Matt Nunes and Dennis Prangle---\n",
"  Current package version: ",ver," (",curdate,") \n\n",
"\n",
"**********************************************\n",
"\n",

"abctools", ver, "loaded\n")
)

}
