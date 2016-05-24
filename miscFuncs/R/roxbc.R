##' roxbc function
##'
##' A function to build and check packages where documentation has been compiled with roxygen. Probably only works in Linux.
##'
##' @param name package name
##' @param checkflags string giving optional check flags to R CMD check, default is --as-cran
##' @return builds and checks the package
##' @export

roxbc <- function(name,checkflags="--as-cran"){ # roxygenize, build and check a package

    roxygenize(name)
    out <- system(paste("R CMD build --compact-vignettes=gs ",name,sep=""),intern=TRUE)
    stuff <- out[length(out)]
    if (length(grep("building",stuff))==0){
        system(paste("R CMD build ",name,sep=""))
        stop("Error in building")
    }
    else{
        buildfn <- substr(stuff,gregexpr("building",stuff)[[1]][1] + 10,nchar(stuff)-1)
        system(paste("R CMD check ",checkflags," ",buildfn,sep=""))
    }
}


##' roxbuild function
##'
##' A function to build  packages where documentation has been compiled with roxygen. Probably only works in Linux.
##'
##' @param name package name
##' @return builds and checks the package
##' @export

roxbuild <- function(name){ # roxygenize, build and check a package

    roxygenize(name)
    out <- system(paste("R CMD build --compact-vignettes=gs ",name,sep=""),intern=TRUE)
    stuff <- out[length(out)]
    if (length(grep("building",stuff))==0){
        system(paste("R CMD build ",name,sep=""))
        stop("Error in building")
    }
}