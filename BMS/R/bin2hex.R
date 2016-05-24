bin2hex <-
function (binvec) 
{
    if (!is.logical(binvec)) {
        if (is.numeric(binvec)) {
            binvec = as.logical(binvec)
        }
        else {
            stop("need to supply a logical vector like c(T,F) or c(1,0)")
        }
    }
    hexobj <- .hexcode.binvec.convert(length(binvec))
    hexcode = hexobj$as.hexcode(binvec)
    if (nchar(hexcode) > (floor((length(binvec) - 1)/4) + 1)) {
        hexcode = substring(hexcode, 2)
    }
    return(hexcode)
}
