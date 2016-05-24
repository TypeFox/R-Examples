c.TextDocument <-
function(..., recursive = FALSE)
{
    args <- list(...)
    x <- args[[1L]]

    if (length(args) == 1L)
        return(x)

    if (!all(unlist(lapply(args, inherits, class(x)))))
        stop("not all arguments are text documents")

    v <- list(content = args,
              meta = CorpusMeta(),
              dmeta = data.frame(row.names = seq_along(args)))
    class(v) <- c("VCorpus", "Corpus")
    v
}

.format_TextDocument <-
function(x, ...)
    c(sprintf("<<%s>>", class(x)[1L]),
      sprintf("Metadata:  %d", length(meta(x))))

PlainTextDocument <-
function(x = character(0),
         author = character(0),
         datetimestamp = as.POSIXlt(Sys.time(), tz = "GMT"),
         description = character(0),
         heading = character(0),
         id = character(0),
         language = character(0),
         origin = character(0),
         ...,
         meta = NULL,
         class = NULL)
{
    p <- list(content = as.character(x),
              meta = TextDocumentMeta(author, datetimestamp, description,
                                      heading, id, language, origin, ...,
                                      meta = meta))
    class(p) <- unique(c(class, "PlainTextDocument", "TextDocument"))
    p
}

as.character.PlainTextDocument <-
function(x, ...)
    content(x)

content.PlainTextDocument <-
function(x)
    x$content

`content<-.PlainTextDocument` <-
function(x, value)
{
    x$content <- as.character(value)
    x
}

format.PlainTextDocument <-
function(x, ...)
    c(.format_TextDocument(x), sprintf("Content:  chars: %d", sum(nchar(x$content))))

meta.PlainTextDocument <-
function(x, tag = NULL, ...)
    if (is.null(tag)) x$meta else x$meta[[tag]]

`meta<-.PlainTextDocument` <-
function(x, tag = NULL, ..., value)
{
    if(is.null(tag))
        x$meta <- value
    else
        x$meta[[tag]] <- value
    x
}

words.PlainTextDocument <-
function(x, ...)
    scan_tokenizer(x)

XMLTextDocument <-
function(x = list(),
         author = character(0),
         datetimestamp = as.POSIXlt(Sys.time(), tz = "GMT"),
         description = character(0),
         heading = character(0),
         id = character(0),
         language = character(0),
         origin = character(0),
         ...,
         meta = NULL)
{
    d <- list(content = x,
              meta = TextDocumentMeta(author, datetimestamp, description,
                                      heading, id, language, origin, ...,
                                      meta = meta))
    class(d) <- c("XMLTextDocument", "TextDocument")
    d
}

as.character.XMLTextDocument <-
function(x, ...)
    as.character(unlist(content(x), use.names = FALSE))

content.XMLTextDocument <-
function(x)
    x$content

`content<-.XMLTextDocument` <-
function(x, value)
{
    x$content <- value
    x
}

format.XMLTextDocument <- .format_TextDocument

meta.XMLTextDocument <- meta.PlainTextDocument

`meta<-.XMLTextDocument` <- `meta<-.PlainTextDocument`
