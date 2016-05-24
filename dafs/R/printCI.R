printCI = function(x, fmt){
    if(length(x)!=2)
        stop("x must have be a vector of length 2")
    strFmt = paste('(',fmt,', ',fmt,')',sep='')
    strResult = sprintf(strFmt, x[1], x[2])
    return(strResult)
}
