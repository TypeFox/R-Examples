`residuals.fm` <- function (object, ...) 
{
    if (class(object)[1] == "fm"|class(object)[1] == "ftsm"){
        return(structure(list(x = object$x1, y = object$y1, z = t(object$residuals$y), xname = object$y$xname, yname = object$y$yname, call = match.call()), 
               class = "fmres"))
    }
    else {
         warning("object is neither a functional time series model nor a functional model.")
    }
}
