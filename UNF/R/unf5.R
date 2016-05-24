unf5 <- 
function(x, 
         digits = 7L, 
         characters = 128L, 
         truncation = 128L,
         raw_as_character = TRUE,
         factor_as_character = TRUE,
         nonfinites_as_missing = FALSE, 
         empty_character_as_missing = FALSE,
         dvn_zero = FALSE,
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
    if(is.complex(x)){
        # COMPLEX numbers: treat as character?
        x <- as.character(x)
        warning("Complex vector converted to character")
    }
    if(is.raw(x)){
        if(raw_as_character) # DVN ingests raw as character
            x <- as.character(x)
        # BIT: Normalize bit fields by converting to big-endian form, truncating all leading empty bits, aligning to a byte boundary by re-padding with leading zero bits, and base64 encoding to form a character string representation.
        else {
            char <- sapply(x, function(i){
                r <- raw()
                as.character(writeBin(i, r, endian='big'))
            })
            warning('UNF is untested on raw vectors')
        }
    }
    if(inherits(x, 'Date')){
        # DATE:
        # Normalize time and date, based on a single, unambiguous representation selected from the many described in the ISO 8601 standard.
        # Convert calendar dates to a character string of the form YYYY-MM-DD. Partial dates in the form YYYY or YYYY-MM are permitted.
        # https://redmine.hmdc.harvard.edu/issues/2997
        if(!date_format %in% c('%Y-%m-%d', '%Y-%m', '%Y', '%F'))
            stop("'date_format' must be '%Y-%m-%d', '%Y-%m', '%Y', or '%F'")
        char <- format(x, fmt = date_format)
    } else if(inherits(x, 'POSIXt')){
        # DATE-TIME: Time representation is based on the ISO 8601 extended format, hh:mm:ss.fffff. When .fffff represents fractions of a second, it must contain no trailing (non-significant) zeroes, and is omitted if valued at zero. Other fractional representations, such as fractional minutes and hours, are not permitted. If the time zone of the observation is known, convert the time value to the UTC time zone and append a "Z" to the time representation.
        if(inherits(x, 'POSIXlt'))
            x <- as.POSIXct(x)
        char <- paste0(format(x, "%Y-%m-%dT%H:%M:", timezone), 
                       gsub("\\.?0+$","",format(x, paste0("%OS",decimal_seconds), timezone)), 
                       ifelse(timezone=="UTC", "Z", ""))
    } else if(is.character(x)){
        # CHARACTER
        char <- as.character(x)
        if(empty_character_as_missing)
            char <- ifelse(x=='',NA,char)
    } else if(is.numeric(x)){
        # NUMERICS: round to nearest, ties to even (use `signif` or `signifz`)
        char <- .expform(signif(x, digits), digits-1)
        if(dvn_zero)
            char <- ifelse(x==0, '+0.e-6', char) # https://redmine.hmdc.harvard.edu/issues/3085
    } else if(is.logical(x)){
        # LOGICAL: normalize boolean to 0, 1, or missing, then treat as numeric
        char <- .expform(as.integer(x), digits-1)
        if(dvn_zero)
            char <- ifelse(x, char, '+0.e-6') # https://redmine.hmdc.harvard.edu/issues/3085
    }
    
    # replace non-finite and missing values with NA
    # https://redmine.hmdc.harvard.edu/issues/2867
    # https://redmine.hmdc.harvard.edu/issues/2960
    out <- .nonfinite(x, char, nonfinites_as_missing, encoding="UTF-8", characters = characters)
    
    hash <- digest(out, algo='sha256', serialize=FALSE, raw=TRUE)
    long <- base64encode(hash)
    short <- base64encode(hash[1:(truncation/8L)]) # truncated UNF
    
    formatted <- paste0('UNF5:',
        if((digits != 7) | (characters != 128)) {
            paste0(paste(digits, characters, sep=','), ':', as.character(short))
        } else {
            as.character(short)
        })
    
    out <- list(unf = as.character(short),
                hash = hash,
                unflong = as.character(long),
                formatted = formatted)
    class(out) <- c('UNF')
    attr(out, 'version') <- 5
    attr(out, 'digits') <- digits
    attr(out, 'characters') <- characters
    attr(out, 'truncation') <- truncation
    return(out)
}
