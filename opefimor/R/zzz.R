.onAttach <- function(libname, pkgname) {
    Pver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
    fields="Version")
    packageStartupMessage(paste(pkgname, Pver))
    packageStartupMessage("Companion package to the book")
    packageStartupMessage(sQuote('Option Pricing and Estimation of Financial Models in R'))
    packageStartupMessage("Wiley, Chichester, S.M. Iacus (2011).")
    packageStartupMessage(paste("For more informations type ",
    sQuote("vignette(\"opefimor\")"),
    "or ?opefimor"))
}
