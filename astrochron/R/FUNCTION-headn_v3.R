### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### headn: head function with column numbers indicated - (SRM: Oct. 31, 2012;
###                                                        May 20, 2013; July 26, 2013)
###
###########################################################################

headn <- function (dat)
{

cat("\n Here are the column numbers associated with each variable in your data frame:\n\n") 

dat <- data.frame(dat)

ncols=length(dat)
header <- head(dat,n=1L)
header[1,] <- 1:ncols

return(header)

### END function headn
}
