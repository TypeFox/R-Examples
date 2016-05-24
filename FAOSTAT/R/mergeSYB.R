##' Function for merging data from different source.
##'
##' This function searches for supported country system and translate
##' the data to allow for join.
##'
##' The names of the data to be merged has to be the same as the
##' FAOcountryProfile code name.
##'
##' @param x data frames, or objects to be coerced to one.
##' @param y data frames, or objects to be coerced to one.
##' @param outCode The country code system to be used to join the
##' different sources.
##' @param all Same as the merge function, defaulted to an outer join.
##' @param ... Arguments to be passed on to the merge function.
##'
##' @export

mergeSYB = function(x, y, outCode = "FAOST_CODE", all = TRUE, ...){
    ## Translate code in x to outCode
    fromCodex = intersect(colnames(x), colnames(FAOcountryProfile))
    if(length(fromCodex) == 0){
        stop("Unknown country code in x, check FAOcountryProfile for supported country code")
    } else if(length(fromCodex) > 1){
        warning(paste("More than one country code system found in x, ",
                      fromCodex[1], " is used.", sep = ""))
        fromCodex = fromCodex[1]
    }
    dfx = translateCountryCode(data = x, from = fromCodex, to = outCode)
    ## Translate code in y to outCode
    fromCodey = intersect(colnames(y), colnames(FAOcountryProfile))
    if(length(fromCodey) == 0){
        stop("Unknown country code in y, check FAOcountryProfile for supported country code")
    } else if(length(fromCodey) > 1){
        warning(paste("More than one country code system found in y, ",
                      fromCodey[1], " is used.", sep = ""))
        fromCodey = fromCodey[1]
    }
    dfy = translateCountryCode(data = y, from = fromCodey, to = outCode)
    ## Merge the translated data
    merge(x = dfx, y = dfy, all = TRUE, ...)
}

utils::globalVariables(names = c("FAOcountryProfile"))