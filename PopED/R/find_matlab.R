## Function written to match MATLAB function
## Author: Andrew Hooker

find_matlab <- function (x) 
{
    expr <- if (is.logical(x)) {
        x
    }
    else {
        x != 0
    }
    return(which(expr))
}

