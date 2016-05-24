`BiodiversityRGUI` <-
function()
{
#    if (! require(vegan)) {stop("Please install the vegan package")}
    options(Rcmdr=list(etc=file.path(path.package(package="BiodiversityR"),
        "etc"), sort.names=FALSE))
    if ("Rcmdr" %in% .packages()) {
        stop("R commander should not have been loaded yet")
    }else{
        Commander()
    }
}

