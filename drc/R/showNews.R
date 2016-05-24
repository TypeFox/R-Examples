"showNews" <- function(pkgname, filename = c("NEWS", "CHANGES"))
{
    filename <- match.arg(filename)

    file.show(paste(.libPaths(), pkgname, filename, sep = "/"), 
    title = paste("Package information for", pkgname))
}
## drc:::showNews("drc")