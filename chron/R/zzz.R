.onLoad <-
function(libname, pkgname) {

    ## The following controls the behavior when faced w. 2-digit years.
    ##
    ## To have 2-digit years actually refer to the first century
    ##    options(chron.year.abb = FALSE)
    ##
    ## To flag all 2-digit years as error:
    ##    options(chron.year.abb = TRUE,
    ##            chron.year.expand = "year.strict")
    ##
    ## To allow 2-digit year abbreviations and guess(?) actual year:
    ##    options(chron.year.abb = TRUE,
    ##            chron.year.expand = "year.expand")

    if(is.null(getOption("chron.year.abb")))
        options(chron.year.abb = TRUE)
    if(is.null(getOption("chron.year.expand")))
        options(chron.year.expand = "year.expand")
}
