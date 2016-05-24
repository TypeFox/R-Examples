#############################################################################################
## File: PKshow.R
## Author: Xiaoyong Sun
## Date: 10/12/2009
## Goal: PKoutput
## Notes:
##      -
#############################################################################################

PKshow <- function(nonmemObj=NULL, table.Rowv=FALSE, table.Colv=FALSE)
{

    PKoutput(nonmemObj, table.Rowv, table.Colv)
    PKcode()
    
    file.location <- paste("file://", getwd(),"/PKindex.html", sep="")
    browseURL(file.location, browser = getOption("browser"))
}

