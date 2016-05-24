#'Hooks for Namespace Events
#'
#'Functions to be called when loaded, attached, detached or unloaded
#'@param libname a character string giving the library directory where
#'@param pkgname a character string giving the name of the package.
.onAttach<-function(libname,pkgname){
    packageStartupMessage("Welcome to package ztable ver 0.1.5")
}

.onLoad<-function(libname,pkgname){
    options(ztable.include.rownames=TRUE)
    options(ztable.include.colnames=TRUE)
    options(ztable.type="viewer")
    options(ztable.color="black")
    options(ztable.show.heading=TRUE)
    options(ztable.show.footer=TRUE)
    options(ztable.caption.placement="top")
    options(ztable.caption.position="c")
    options(ztable.caption.bold=FALSE)
    options(ztable.booktabs=FALSE)
    options(ztable.zebra=NULL)
    options(ztable.zebra.color=NULL)
    options(ztable.zebra.type=1)
    options(ztable.zebra.colnames=FALSE)
    options(ztable.zebra.rownames=TRUE)
    options(ztable.colnames.bold=FALSE)
    invisible()
}
