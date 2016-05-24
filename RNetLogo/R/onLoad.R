.rnetlogo <- new.env()
.rnetlogo[["savedworkingdir"]] <- c()
.rnetlogo[["libname"]] <- ""
.rnetlogo[["pkgname"]] <- ""
.rnetlogo[["nlversion"]] <- 0
.rnetlogo[["nl3d"]] <- -1
.rnetlogo[["guiobj"]] <- NULL
.rnetlogo[["objects"]] <- c()

.onLoad <-
function(libname, pkgname)
{
    .rnetlogo$pkgname <- pkgname 	
    .rnetlogo$libname <- libname
    .rnetlogo$startedGUI <- FALSE
}
    