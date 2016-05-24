maybeQuote <- function (x) {
    if (length(x) == 0) 
        return(character())
    if(is.numeric(x))return(x)
    #if (l10n_info()$"UTF-8") 
     #   paste("?", x, "?", sep = "")
    #else 
	paste("'", x, "'", sep = "")
}

