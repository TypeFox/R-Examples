`print.fds` <- function (x, digits = NULL, ...) 
{
    if (class(x)[1] == "fds"){
        cat("Functional data set")
        cat(paste("\n y:", x$yname))
        cat(paste("\n x:", x$xname, "\n"))
    }
    else {
         stop("object is not a functional data set.")
    }
}

