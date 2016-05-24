## last modified 23 Oct 2002 by J. Du

is.mixdata <- function (x) 
{
    flag <- TRUE
    if (sum(!is.na(match(class(x), "mixdata"))) == 0) 
        flag <- FALSE
    flag
}
