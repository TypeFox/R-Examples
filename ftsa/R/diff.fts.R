`diff.fts` <- function (x, lag = 1, differences = 1, ...) 
{
    if (class(x)[1] == "fts"|class(x)[1] == "sfts"){
        x$y <- t(diff(t(x$y), lag, differences, ...))
        return(x)
    }
    else{
        stop("object is not a functional time series or a sliced functional time series")
    }
}

