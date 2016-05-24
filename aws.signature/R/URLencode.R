URLencode <- function(URL, reserved = FALSE) { 
    OK <- paste0("[^", if (!reserved) 
        "][!$&'()*+,;=:/?@#", "ABCDEFGHIJKLMNOPQRSTUVWXYZ", 
        "abcdefghijklmnopqrstuvwxyz0123456789._~-", "]") 
    x <- strsplit(URL, "")[[1L]] 
    z <- grep(OK, x) 
    if (length(z)) { 
        y <- sapply(x[z], function(x) paste0("%", toupper(as.character(charToRaw(x))), collapse = "")) 
        x[z] <- y 
    } 
    paste(x, collapse = "") 
} 
