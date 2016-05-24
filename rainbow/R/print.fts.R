`print.fts` <- function (x, digits = NULL, ...) 
{
    if (class(x)[1] == "fts"){
        cat("Functional time series")
        cat(paste("\n y:", x$yname))
        cat(paste("\n x:", x$xname, "\n"))
    }
    else {
         stop("object is not a functional time series.")
    }
}

