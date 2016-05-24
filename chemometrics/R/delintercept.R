delintercept <- function (mm) 
{
    saveattr <- attributes(mm)
    intercept <- which(saveattr$assign == 0)
    if (!length(intercept)) 
        return(mm)
    mm <- mm[, -intercept, drop = FALSE]
    saveattr$dim <- dim(mm)
    saveattr$dimnames <- dimnames(mm)
    saveattr$assign <- saveattr$assign[-intercept]
    attributes(mm) <- saveattr
    mm
}

