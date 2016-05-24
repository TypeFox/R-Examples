# Author: Ingo Feinerer

TextDocumentMeta <-
function(author, datetimestamp, description, heading, id, language, origin, ...,
         meta = NULL)
{
    if (is.null(meta))
        meta <- list(author = author, datetimestamp = datetimestamp,
                     description = description, heading = heading, id = id,
                     language = language, origin = origin, ...)

    stopifnot(is.list(meta))
    if (!is.null(meta$author) && !inherits(meta$author, "person"))
        meta$author <- as.character(meta$author)
    if (!is.null(meta$datetimestamp) && !inherits(meta$datetimestamp, "POSIXt"))
        meta$datetimestamp <- as.character(meta$datetimestamp)
    if (!is.null(meta$description))
       meta$description <- as.character(meta$description)
    if (!is.null(meta$heading))
       meta$heading <- as.character(meta$heading)
    if (!is.null(meta$id))
       meta$id <- as.character(meta$id)
    if (!is.null(meta$language))
       meta$language <- as.character(meta$language)
    if (!is.null(meta$origin))
       meta$origin <- as.character(meta$origin)

    class(meta) <- "TextDocumentMeta"
    meta
}

print.TextDocumentMeta <-
function(x, ...)
{
    cat(sprintf("  %s: %s",
                format(names(x), justify = "left"),
                sapply(x, as.character)),
        sep = "\n")
    invisible(x)
}

CorpusMeta <-
function(..., meta = NULL)
{
    if (is.null(meta))
        meta <- list(...)

    stopifnot(is.list(meta))

    class(meta) <- "CorpusMeta"
    meta
}

meta.VCorpus <- meta.PCorpus <-
function(x, tag = NULL, type = c("indexed", "corpus", "local"), ...)
{
    if (!is.null(tag) && missing(type)) {
        type <- if (tag %in% colnames(x$dmeta)) "indexed"
        else if (tag %in% names(x$meta)) "corpus"
        else "local"
    }
    type <- match.arg(type)
    if (identical(type, "indexed"))
        if (is.null(tag)) x$dmeta else x$dmeta[tag]
    else if (identical(type, "corpus"))
        if (is.null(tag)) x$meta else x$meta[[tag]]
    else if (identical(type, "local"))
        lapply(x, meta, tag)
    else
        stop("invalid type")
}

`meta<-.VCorpus` <- `meta<-.PCorpus` <-
function(x, tag, type = c("indexed", "corpus", "local"), ..., value)
{
    type <- match.arg(type)
    if (identical(type, "indexed"))
        x$dmeta[, tag] <- value
    else if (type == "corpus")
        x$meta[[tag]] <- value
    else if (identical(type, "local")) {
        for (i in seq_along(x))
            meta(x[[i]], tag) <- value[i]
    } else
        stop("invalid type")
    x
}

# Simple Dublin Core to tm metadata mapping
# http://en.wikipedia.org/wiki/Dublin_core#Simple_Dublin_Core
Dublin_Core_tm_map <-
list("contributor" = "contributor",
     "coverage" = "coverage",
     "creator" = "author",
     "date" = "datetimestamp",
     "description" = "description",
     "format" = "format",
     "identifier" = "id",
     "language" = "language",
     "publisher" = "publisher",
     "relation" = "relation",
     "rights" = "rights",
     "source" = "source", # or better "origin"?
     "subject" = "subject",
     "title" = "heading",
     "type" = "type"
     )

DublinCore <-
function(x, tag = NULL)
{
    tmm <- unlist(Dublin_Core_tm_map, use.names = FALSE)
    dcm <- names(Dublin_Core_tm_map)

    if (is.null(tag)) {
        m <- lapply(tmm, function(t) meta(x, t))
        names(m) <- dcm
        class(m) <- "TextDocumentMeta"
        m
    } else
        meta(x, tmm[charmatch(tolower(tag), dcm)])
}

`DublinCore<-` <-
function(x, tag, value)
{
    tmm <- unlist(Dublin_Core_tm_map, use.names = FALSE)
    dcm <- names(Dublin_Core_tm_map)

    meta(x, tmm[charmatch(tolower(tag), dcm)]) <- value
    x
}
