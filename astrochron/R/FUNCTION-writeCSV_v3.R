### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### writeCSV function - (SRM: Feb. 15, 2012; May 20, 2013; May 1, 2014)
###
### output a CSV file
###########################################################################

writeCSV <- function (filename,output)
{

 cat("\n----- WRITING CSV FILE -----\n")
 write.table(file=filename, output, sep=",",row.names=FALSE)

### END function writefile
}
