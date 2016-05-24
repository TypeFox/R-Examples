

.onAttach <- function(libname, pkgname) {
    RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
	packageStartupMessage(" ")
    packageStartupMessage(paste(pkgname, RFver))
	packageStartupMessage(" ")
    packageStartupMessage("Copyright (C) 2006 - 2012, Alejandro Jara")
    packageStartupMessage("Department of Statistics")
    packageStartupMessage("P.U. Catolica de Chile")
    packageStartupMessage(" ")
    packageStartupMessage("Support provided by Fondecyt")
    packageStartupMessage("11100144 grant.")
    packageStartupMessage(" ")
}


.onUnload <- function(libpath) {
    library.dynam.unload("DPpackage", libpath)
}


"print.anovaPsCP"<-
function (x, digits = max(getOption("digits") - 2, 3), signif.stars = getOption("show.signif.stars"), 
    ...) 
{
    if (!is.null(heading <- attr(x, "heading"))) 
        cat(heading, sep = "\n")

    nc <- dim(x)[2]
    if (is.null(cn <- colnames(x))) 
        stop("'anova' object must have colnames")

    has.P <- substr(cn[nc], 1, 3) == "PsC"

    printCoefmat(x, digits = digits, signif.stars = signif.stars, 
        has.Pvalue = has.P, P.values = has.P,eps.Pvalue=0.01)
    invisible(x)
}

