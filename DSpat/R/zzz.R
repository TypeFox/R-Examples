print.DSpat.version <- function()
{ library(help=DSpat)$info[[1]] -> version
	version <- version[pmatch("Version",version)]
	if(!is.null(version)){
		um <- strsplit(version," ")[[1]]
		version <- um[nchar(um)>0][2]
	}	
	hello <- paste("This is DSpat ",version,"\n",sep="")
	packageStartupMessage(hello)
}

.onAttach <- function(...)  
{
	print.DSpat.version()
}

