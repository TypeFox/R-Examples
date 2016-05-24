
read.pattern <- function(file, pattern = "[^[:space:]]+", perl = FALSE, text, sep = "\01", fileEncoding = "", ...) {

    if (missing(file) && !missing(text)) {
        file <- textConnection(text)
        on.exit(close(file))
    }
    if (is.character(file)) {
        file <- if (nzchar(fileEncoding)) 
            file(file, "rt", encoding = fileEncoding)
        else file(file, "rt")
        on.exit(close(file))
    }
    if (!inherits(file, "connection")) 
        stop("'file' must be a character string or connection")
    if (!isOpen(file, "rt")) {
        open(file, "rt")
        on.exit(close(file))
    }
    text <- readLines(file)
    tmp <- sapply(strapplyc(text, pattern), paste, collapse = "\01")
    read.table(text = tmp, sep = sep, ...)

}

