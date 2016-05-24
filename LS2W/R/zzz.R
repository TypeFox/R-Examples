## this is the equivalent of .First.lib (for packages with a namespace)
## note that it is general, so the name of the package etc doesn't need
## to be specified (just say whatever you want in the prints though).

.onAttach <-function(lib,pkg)
{
ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
     ver <- as.character(ver)	
curdate <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Date")
    curdate <- as.character(curdate)	

# Welcome message (MAN):

packageStartupMessage(paste(
"\n",
"***********************************************************\n",
"LS2W: a package for 2D Locally Stationary Wavelet processes\n\n",
"      --- Written by Idris Eckley and Guy Nason ---\n",
"  --- Contributions from Sarah Taylor and Matt Nunes ---\n",
"   Current package version: ",ver," (",curdate,") \n\n",
"\n",
"***********************************************************\n",
"\n")
)

}

# mode argument below not really necessary:
if(!exists("DWEnv",mode="environment")){
    DWEnv<-new.env()
} 
