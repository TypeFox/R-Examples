WordListDocument <-
function(con, encoding = "unknown", meta = list())
{
    words <- readLines(con, encoding = encoding, warn = FALSE)
    doc <- list(content = words, meta = meta)
    class(doc) <- c("WordListDocument", "TextDocument")
    doc
}

format.WordListDocument <-
function(x, ...)    
    c(.format_TextDocument(x),
      sprintf("Content:  words: %d", length(x$content)))
    
## print.WordListDocument <-
## function(x, ...)
## {
##     writeLines(sprintf("<<WordListDocument (words: %d)>>",
##                        length(x$content)))
##     invisible(x)
## }

content.WordListDocument <-
function(x)
    x$content

## meta.WordListDocument <-
## function(x, tag = NULL, ...)
##     if(is.null(tag)) x$meta else x$meta[[tag]]

## `meta<-.WordListDocument` <-
## function(x, tag = NULL, ..., value)
## {
##     if(is.null(tag))
##         x$meta <- value
##     else
##         x$meta[[tag]] <- value
##     x
## }

as.character.WordListDocument <-
words.WordListDocument <-
function(x, ...)
    x$content
