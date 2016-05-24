.First.lib <-
    function (libname, pkgname)
{
    ## echo output to screen
    cat("## Statistical Methods in Forensic Genetics \n")
    cat("## The work was supported by the project 1M06014 of the Ministry of 
        Education, Youth and Sports of the Czech Republic. \n")
    ## assuming you would need the package genetics, combinat for your package
    require(genetics, quietly=TRUE)
    require(combinat, quietly=TRUE)
}
