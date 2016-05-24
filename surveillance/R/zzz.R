#######################################
### Hook functions for package start-up
#######################################

gpclibCheck <- function (fatal = TRUE)
{
    gpclibOK <- surveillance.options("gpclib")
    if (!gpclibOK && fatal) {
        message("Note: The gpclib license is accepted by ",
                sQuote("surveillance.options(gpclib=TRUE)"), ".")
        stop("acceptance of the gpclib license is required")
    }
    gpclibOK
}

.onLoad <- function (libname, pkgname)
{
    ## initialize options
    reset.surveillance.options()
}

.onAttach <- function (libname, pkgname)
{
    ## Startup message
    VERSION <- packageVersion(pkgname, lib.loc=libname)
    packageStartupMessage("This is ", pkgname, " ", VERSION, ". ",
                          "For overview type ",
                          sQuote(paste0("help(", pkgname, ")")), ".")

    ## decide if we should run all examples (some take a few seconds)
    allExamples <- if (interactive()) {
        TRUE
    } else { # R CMD check
        ## only do all examples if a specific environment variable is set
        ## (to any value different from "")
        nzchar(Sys.getenv("_R_SURVEILLANCE_ALL_EXAMPLES_"))
        ## CAVE: testing for _R_CHECK_TIMINGS_ as in surveillance < 1.9-1
        ## won't necessarily skip long examples for daily checks on CRAN (see
        ## https://stat.ethz.ch/pipermail/r-devel/2012-September/064812.html
        ## ). For instance, the daily Windows checks run without timings.
    }
    surveillance.options(allExamples = allExamples)
}



###########################
### Little helper functions
###########################


### checking if x is scalar, i.e. a numeric vector of length 1.

isScalar <- function (x) {
    length(x) == 1L && is.vector(x, mode = "numeric")
}


### _c_onditional lapply, which only uses lapply() if X really is a list object
### and otherwise applies FUN to X. The result is always a list (of length 1 in
### the latter case). Used for neOffset in hhh4 models.

clapply <- function (X, FUN, ...)
{
    if (is.list(X)) lapply(X, FUN, ...) else list(FUN(X, ...))
}


### pretty p-value formatting

formatPval <- function (pv, eps = 1e-4, scientific = FALSE, ...)
{
    format1 <- function (p)
        format.pval(p, digits = if (p < 10*eps) 1 else 2, eps = eps,
                    nsmall = 2, scientific = scientific, ...)
    vapply(X = pv, FUN = format1, FUN.VALUE = "", USE.NAMES = TRUE)
}


### determines multiplicities in a matrix (or data frame)
### and returns unique rows with appended column of counts

countunique <- function (x) unique(cbind(x, COUNT = multiplicity(x)))


### generate a color vector (via the colorspace package)

hcl.colors <- function (ncolors=100, use.color=TRUE)
{
    GYR <- if (requireNamespace("colorspace", quietly=TRUE)) {
        ## the Zeil-ice colors 
        colorspace::heat_hcl(ncolors, h=c(0,120),
                             c=if (use.color) c(90,30) else c(0,0),
                             l=c(50,90), power=c(0.75, 1.2))
    } else {
        if (use.color) heat.colors(ncolors) else grey.colors(ncolors)
    }
    
    return(rev(GYR))
}
