`print.sfts` <- function (x, ...) 
{
    if (class(x)[1] == "sfts"){
        cat("Sliced functional time series")
        cat(paste("\n y:", x$yname))
        cat(paste("\n x:", x$xname, "\n"))
    }
    else {
         stop("object is not a sliced functional time series.")
    }
}

