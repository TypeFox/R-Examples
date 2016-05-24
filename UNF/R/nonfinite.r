.nonfinite <- 
function(x, char, 
         nonfinites_as_missing = FALSE, 
         encoding = "",
         characters = 128L){
    if(nonfinites_as_missing){
        char <- ifelse(!is.finite(x), NA, char)
    } else {
        char <- ifelse((!is.finite(x) & !is.character(x)), tolower(as.character(x)), char)
        char <- ifelse(char=='inf', '+inf', char)
        char <- ifelse(is.nan(x), '+nan', char)
    }
    char <- ifelse(is.na(x) & !is.nan(x), NA, char)
    unicode <- iconv(substring(char,1,characters), to = encoding, toRaw=TRUE)
    out <- unlist(lapply(unicode, function(i) {
        if(is.null(i)) {
            # append three null bits
            intToBits(0)[1:3] 
        } else {
            # append EOL and null bit
            c(i,charToRaw("\n"),intToBits(0)[1])
        }
    }))
    return(out)
}
