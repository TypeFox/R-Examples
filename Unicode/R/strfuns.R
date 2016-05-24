Unicode_alphabetic_tokenizer <-
function(x)
{
    x <- .str_to_u_char(paste(x, collapse = "\n"))
    x[u_char_property(x, "Alphabetic") == "No"] <- 32L
    unlist(strsplit(intToUtf8(x), " +"))
}

.str_to_u_char <-
function(x)
{
    x <- as.character(x)
    if(length(x) > 1L) {
        if(any(nchar(x) != 1L))
            stop("Invalid 'x'.")
        x <- paste(x, collapse = "")
    }
    if(Encoding(x) == "latin1")
        x <- iconv(x, from = "latin1", to = "UTF-8")
    x <- if(is.na(x))
        NA_integer_
    else
        utf8ToInt(x)
    as.u_char(x)
}
