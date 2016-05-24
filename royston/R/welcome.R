.onAttach = function(libname, pkgname) {
    msg = "The Royston's H test has been moved to 'MVN' package. 'royston' package will no longer be supported. Please use 'MVN' package for further analysis."
    msg = strwrap(msg, exdent=4, indent=4)
    packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}
