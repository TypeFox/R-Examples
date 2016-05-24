oauth.encode1 <-
function(x) {
    encode <- function(x) paste0("%", toupper(as.character(charToRaw(x))))
    
    x <- as.character(x)
    chars <- strsplit(x, "")[[1]]
    ok <- !str_detect(chars, "[^A-Za-z0-9_.~-]")
    
    if (all(ok)) return(x)
    
    chars[!ok] <- unlist(lapply(chars[!ok], encode))
    paste0(chars, collapse = "")
}
