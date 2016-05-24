PDF_fonts <-
function(file)
{
    x <- PDF_doc_from_file(file)

    y <- .Call("PDF_fonts", x$doc)

    errors <- attr(y, "errors")
    
    n <- length(y)
    if(n > 0L) {
        y <- unlist(rev(y), recursive = FALSE, use.names = FALSE)
        k <- length(y) / n
        y <- split(y, rep.int(seq_len(k), n))

        y <- as.data.frame(lapply(y[-1L], unlist),
                           stringsAsFactors = FALSE)
    } else {
        y <- data.frame(character(), character(), character(),
                        logical(), logical(),
                        stringsAsFactors = FALSE)
    }

    names(y) <- c("name", "type", "file", "emb", "sub")
    class(y) <- c("PDF_fonts", class(y))
    attr(y, "errors") <- errors

    y
}

print.PDF_fonts <-
function(x, ...)
{
    NextMethod("print", x, ...)
    if(length(errors <- attr(x, "errors")))
        writeLines(c("\nPDF problems:", sprintf("  %s", errors)))
    invisible(x)
}    
