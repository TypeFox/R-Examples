sessionInfoX <- function(pkgs=NULL, list.libP = FALSE, extraR.env = TRUE) {
    ## return an object; then print() via method
    if(!is.null(pkgs)) stopifnot(is.character(pkgs), length(pkgs) > 0)
    lP <- .libPaths() # *is* normalized in the sense of normalizePath()
    nRL <- normalizePath(RLIBS <- strsplit(Sys.getenv("R_LIBS"), ":")[[1]])
    structure(class = "sessionInfoX",
        list(sInfo = sessionInfo(),
             pkgDescr = if(!is.null(pkgs)) sapply(pkgs, packageDescription, simplify=FALSE),
             libPath = lP, .Library = .Library, RLIBS = RLIBS, n.RLIBS = nRL,
             list.libP = if(list.libP) sapply(lP, list.files, simplify=FALSE),
             R.env = Sys.getenv(c("R_ENVIRON", "R_PROFILE", "R_CHECK_ENVIRON")),
             xR.env = if(extraR.env) local({
                 ss <- Sys.getenv()
                 ss[grepl("^_?R_", names(ss))]
             })))
}

print.sessionInfoX <- function(x, ...) {
    cat("Extended  sessionInfo():",
        "-----------------------", sep="\n")# does add a final '\n'
    if(!is.null(pkgD <- x$pkgDescr)) {
        cat("specific  packageDescription()s:\n")
        print(pkgD, ...)
        cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    }
    cat("R_LIBS:\n")
    cbind(x$RLIBS)
    xtr.lp <- setdiff(x$libPath,
                      union(normalizePath(x$.Library), x$n.RLIBS))
    if(length(xtr.lp)) {
        cat("libPath contents in addition to RLIBS and .Library:\n")
        print(xtr.lp)
    } else
        cat("libPath contains RLIBS and .Library (normalized)\n")
    if(length(xx <- setdiff(x$n.RLIBS, x$libPath))) { ## typically empty
        cat("** RLIBS has entries not in .libPaths():\n")
        print(xx)
    }
    cat("Main R env. variables (for more, inspect the 'xR.env' component):\n")
    print(cbind(x$R.env), ...)
    cat("---------------- standard sessionInfo():\n")
    print(x$sInf, ...)
    invisible(x)
}
