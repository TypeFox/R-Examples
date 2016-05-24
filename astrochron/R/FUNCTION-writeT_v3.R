### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### writeT function - (SRM: Mar. 26, 2012; May 20, 2013; May 1, 2014)
###
### output a tab-delimited text file
###########################################################################

writeT <- function (filename,output)
{

 cat("\n----- WRITING TAB-DELIMITED FILE -----\n")
 write.table(file=filename, output, sep="\t",row.names=FALSE)

### END function writeT
}
