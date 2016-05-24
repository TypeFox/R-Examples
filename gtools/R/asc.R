asc <- function(char, simplify=TRUE)
    sapply(char, function(x) strtoi(charToRaw(x),16L), simplify=simplify )

chr <- function(ascii) sapply(ascii, function(x) rawToChar(as.raw(x)) )
