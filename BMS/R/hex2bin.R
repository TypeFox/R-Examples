hex2bin <-
function (hexcode) 
{
    if (!is.character(hexcode)) 
        stop("please input a character like '0af34c'")
    hexcode <- paste("0", tolower(hexcode), sep = "")
    hexobj <- .hexcode.binvec.convert(length(hexcode) * 16L)
    return(hexobj$as.binvec(hexcode))
}
