#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* campus *.* lmu *.* de]
## Time-stamp: <[getNumbers.R] by DSB Die 07/07/2009 14:49 (CEST)>
##
## Description:
## Internal helper function for extractSamples, which extracts numbers from (BayesX
## log file) strings.
#####################################################################################


getNumbers <- function(beforeStringsList,  # list with the strings before the numbers
                       stringVector) # look for the strings in this vector
{
    lapply(beforeStringsList,
           function(string)
           as.numeric(sub(pattern=".+[[:blank:]]+([[:digit:]]+)[[:blank:]]*$",
                          replacement="\\1",
                          x=
                          stringVector[grep(pattern=
                                            paste(".*",
                                                  string,
                                                  sep=""),
                                            x=stringVector)])))
}

