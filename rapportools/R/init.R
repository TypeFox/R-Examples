.onLoad <- function(libname, pkgname)
{

    ## use labels if not set before by the rapport package or user
    if (is.null(getOption('rapport.use.labels')))
        options('rapport.use.labels' = TRUE)

}
