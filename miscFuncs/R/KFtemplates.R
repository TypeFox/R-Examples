##' KFtemplates function
##'
##' A function to print KFfit and KFparest templates to the console. See vignette("miscFuncs")
##' for more information
##'
##' @return Tust prints to the console. This can be copied and pasted into a text editor for further manipulation. 
##' @export

KFtemplates <- function(){
    path <- file.path(find.package("miscFuncs", lib.loc=NULL, quiet = TRUE),"inst")
    fitpath <- file.path(path,"KFfit.R")
    parpath <- file.path(path,"KFparest.R")
    
    l <- readLines(fitpath)
    for(i in 1:length(l)){
        cat(l[[i]],"\n")
    }
    l <- readLines(parpath)
    for(i in 1:length(l)){
        cat(l[[i]],"\n")
    }
    cat("\n\n\n\n")
}

