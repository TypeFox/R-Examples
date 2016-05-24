###
### $Id: zzz.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Hooks called as package is loaded, attached, detatched, and unloaded.
###


## Namespace environment for this package
.MatlabNamespaceEnv <- new.env()


##
## Package/Namespace Hooks
##

##-----------------------------------------------------------------------------
.onAttach <- function(libname, pkgname) {
    verbose <- getOption("verbose")
    if (verbose) {
        local({
            libraryPkgName <- function(pkgname, sep="_") {
                unlist(strsplit(pkgname, sep, fixed=TRUE))[1]
            }
            packageDescription <- function(pkgname) {
                fieldnames <- c("Title", "Version")
                metafile <- file.path(libname, pkgname, "DESCRIPTION")
                meta <- as.list(read.dcf(metafile, fieldnames))
                names(meta) <- fieldnames
                return(meta)
            }

            meta <- packageDescription(pkgname)
            msg <- sprintf("%s, version %s",
                           meta$Title, meta$Version)
            packageStartupMessage(msg)
            msg <- sprintf("Type library(help=%s) to see package documentation",
                           libraryPkgName(pkgname))
            packageStartupMessage(msg)
        })
    }
}


##-----------------------------------------------------------------------------
.onLoad <- function(libname, pkgname) {
    ## Create environment for internal usage
    environment(.MatlabNamespaceEnv) <- asNamespace("matlab")

    ## Load internal variables
    assign("savedTime", 0, envir=.MatlabNamespaceEnv)

    ## Allow no changes or additions to environment
    lockEnvironment(.MatlabNamespaceEnv, bindings=TRUE)

    ## Only allow internal vars to change
    unlockBinding("savedTime", .MatlabNamespaceEnv)
}

