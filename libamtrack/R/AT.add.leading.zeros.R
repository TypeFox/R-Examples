AT.add.leading.zeros <- function(x, digits = 5){
# Convert 99 into "00099" etc.
# From the s-mailing list.
# [S] SUMMARY:  Formatting with Leading Zeroes
# June 21, 2001, Mike.Prager@noaa.gov
# The following solution is based on the one suggested by Bill Venables.
# Library: clan
# Created July 9, 2001, Revised: July 9, 2001, Claus E. Andersen
# Modified for R use (penalty on scientific notation), Dec 1, 2010, Steffen Greilich
options(scipen=3)
mx <- paste("1",paste(rep("0",digits),sep="",collapse=""),sep="")
return(substring(as.numeric(mx)+x, 2))
options(scipen=0)
}