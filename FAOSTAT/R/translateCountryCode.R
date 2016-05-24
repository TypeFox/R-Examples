##' A function to translate between different country coding systems
##'
##' The function translate any country code scheme to another if both
##' are in the FAOcountryProfile
##'
##' @param data The data frame
##' @param from The name of the old coding system
##' @param to The name of the new coding system
##' @param oldCode The column name of the old country coding scheme
##' @export
##'


translateCountryCode = function (data, from, to, oldCode)
{
    cat("\nNOTE: Please make sure that the country are matched according to their definition\n\n")
    if (missing(oldCode))
        oldCode = from
    if (from != to) {
        codeTrans = FAOcountryProfile[which(FAOcountryProfile[,
            from] %in% data[, oldCode]), c(from, to)]
        trans.df = merge(x = codeTrans, y = data, by.x = from,
            by.y = oldCode, all.y = TRUE)
        if (any(is.na(trans.df[, to]))) {
            warning(paste("The following entries does not have '",
                to, "' available\n", sep = ""), immediate. = TRUE)
            print(FAOcountryProfile[FAOcountryProfile[, from] %in%
                                    trans.df[is.na(trans.df[, to]), from],
                                    c(from, to, "OFFICIAL_FAO_NAME")])
        }
        trans.df$OFFICIAL_FAO_NAME = NULL
    }
    else {
        trans.df = data
    }
    trans.df
}

utils::globalVariables(names = c("FAOcountryProfile"))

## translateCountryCode = function(data, from, to, oldCode){
##     warning("Please make sure that the country are matched according to their definition")
##     if(missing(oldCode))
##         oldCode = from
##     if(from != to){
##         codeTrans = FAOcountryProfile[which(FAOcountryProfile[, from] %in%
##         data[, oldCode]), c(from, to, "OFFICIAL_FAO_NAME")]
##         trans.df = merge(x = codeTrans, y = data, by.x = from,
##             by.y = oldCode, all.y = TRUE)
##         if(any(is.na(trans.df[, to]))){
##             warning(paste("The following entries does not have '",
##                           to, "' available", sep = ""))
##             print(unique(trans.df[is.na(trans.df[, to]),
##                                   c(from, to, "OFFICIAL_FAO_NAME")]))
##         }
##         trans.df$OFFICIAL_FAO_NAME = NULL
##     } else {
##         trans.df = data
##     }
##     trans.df
## }
