## R/S-Plus compatibility
## $Id: S4R.R 189 2011-10-01 13:16:39Z dirk.eddelbuettel $

## This package was developed as a part of Summer of Code program organized by Google.
## Thanks to David A. James & Saikat DebRoy, the authors of RMySQL package.
## Code from RMySQL package was reused with the permission from the authors.
## Also Thanks to my GSoC mentor Dirk Eddelbuettel for helping me in the development.


usingR <- function(major=0, minor=0){
    if(is.null(version$language))
        return(FALSE)
    if(version$language!="R")
        return(FALSE)
    version$major>=major && version$minor>=minor
}

## constant holding the appropriate error class returned by try()
if(usingR()){
    ErrorClass <- "try-error"
} else {
    ErrorClass <- "Error"
}


