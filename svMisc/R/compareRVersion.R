compareRVersion <- function (version)
{
    ## This is similar to compareVersion, but works for R version comparison
    compareVersion(paste(R.version$major, R.version$minor, sep = "."),
		version)
}
