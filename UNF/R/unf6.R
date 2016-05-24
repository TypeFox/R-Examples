unf6 <-
function(x, 
         digits = 7L, 
         characters = 128L, 
         truncation = 128L,
         raw_as_character = TRUE,
         factor_as_character = TRUE,
         complex_as_character = TRUE,
         nonfinites_as_missing = FALSE, 
         timezone = "",
         date_format = "%Y-%m-%d",
         decimal_seconds = 5,
         ...){
    if(!truncation %in% c(128,192,196,256))
        stop("'truncation' must be in 128, 192, 196, 256")
    if(truncation < characters)
        stop("'truncation' must be greater than or equal to 'characters'")
    if(inherits(x, 'AsIs')){
        x <- as.character(x)
    }
    if(inherits(x, 'ts') | inherits(x, 'zoo') | inherits(x, 'difftime')) {
        x <- as.numeric(x)
    }
    if(is.factor(x)) {
        # FACTOR: treat factor as character and truncate to k
        if(factor_as_character)
            x <- as.character(x)
        else
            x <- as.numeric(x)
    }
    if(is.raw(x)){
        if(raw_as_character) # DVN ingests raw as character
            x <- as.character(x)
        # BIT: Normalize bit fields by converting to big-endian form,
        #      truncating all leading empty bits, 
        #      aligning to a byte boundary by re-padding with leading zero bits, and 
        #      base64 encoding to form a character string representation
        else {
            char <- sapply(x, function(i){
                r <- raw()
                as.character(writeBin(i, r, endian='big'))
            })
            warning('UNF is untested on raw vectors')
        }
    }
    if(is.complex(x) & complex_as_character) {
        x <- as.character(x)
    }
    if(is.complex(x) & !complex_as_character){
        # COMPLEX numbers: format as `A,iB`
        re <- .expform(signif(Re(x), digits), digits-1)
        co <- .expform(signif(Im(x), digits), digits-1)
        char <- paste(re, co, sep=",i")
    } else if(inherits(x, 'Date')){
        # DATE: Dates are converted to character strings in the form "YYYY-MM-DD", but partial dates ("YYYY" and "YYYY-MM") are permitted.
        if(!date_format %in% c('%Y-%m-%d', '%Y-%m', '%Y', '%F'))
            stop("'date_format' must be '%Y-%m-%d', '%Y-%m', '%Y', or '%F'")
        char <- format(x, fmt = date_format)
    } else if(inherits(x, 'POSIXt')){
        # DATE-TIME: Encoded as: `"2014-08-22T16:51:05Z"`.
        if(inherits(x, 'POSIXlt'))
            x <- as.POSIXct(x)
        char <- paste0(format(x, "%Y-%m-%dT%H:%M:", timezone), 
                       gsub("\\.?0+$","",format(x, paste0("%OS",decimal_seconds), timezone)), 
                       ifelse(timezone=="UTC", "Z", ""))
    } else if(is.numeric(x)){
        # NUMERICS: round to nearest, ties to even (use `signif` or `signifz`)
        char <- .expform(signif(x, digits), digits-1)
    } else if(is.logical(x)){
        # LOGICAL: normalize boolean to 0, 1, or missing, then treat as numeric
        char <- .expform(as.integer(x), digits-1)
    } else {
        char <- as.character(x)
    }
    
    # deal with non-finite and missing values; convert to raw
    out <- .nonfinite(x, char, nonfinites_as_missing, encoding = "UTF-8", characters = characters)
    
    hash <- digest(out, algo='sha256', serialize=FALSE, raw=TRUE)
    long <- base64encode(hash)
    short <- base64encode(hash[1:(truncation/8L)]) # truncated UNF
    
    # format printable UNF
    header <- paste(if(digits != 7) paste0("N", digits) else NULL,
                    if(characters != 128) paste0("X", characters) else NULL,
                    if(truncation != 128) paste0("H", truncation) else NULL,
                    sep = ",", collapse="")
    header <- ifelse(length(header), gsub("^[[:punct:]]+", "", header), "")
    header <- ifelse(length(header), gsub("[[:punct:]]+$", "", header), "")
    header <- ifelse(length(header), gsub("[[:punct:]]{2}", ",", header), "")
    formatted <- paste0('UNF6:', ifelse(header == "", 
                                        as.character(short), 
                                        paste0(header,':', as.character(short))))
    out <- list(unf = as.character(short),
                hash = hash,
                unflong = as.character(long),
                formatted = formatted)
    class(out) <- c('UNF')
    attr(out, 'version') <- 6
    attr(out, 'digits') <- digits
    attr(out, 'characters') <- characters
    attr(out, 'truncation') <- truncation
    return(out)
}
