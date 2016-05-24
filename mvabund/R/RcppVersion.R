RcppVersion <- function() {

    version <- installed.packages()["Rcpp","Version"]

    ## Print version of this package
    cat('\nRcpp package Version ', version, '\n')

    ## Print version of Rcpp/RcppTemplate used to build this package
    licenseFile <- system.file(".", "LICENSE-Rcpp.txt", package="Rcpp")
    con <- file(licenseFile, "r")
    writeLines(readLines(con)[1:4])
    close(con)
    invisible(NULL)
}
