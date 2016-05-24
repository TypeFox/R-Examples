AnnotatedPlainTextDocument <-
function(s, annotations, meta = list())
{
    s <- as.String(s)

    ## Be nice.
    if(is.Annotation(annotations))
        annotations <- list(annotations)
    else {
        annotations <- as.list(annotations)
        if(!length(annotations) ||
           !all(sapply(annotations, is.Annotation)))
            stop("argument 'annotations' must give a positive number of Annotation objects")
    }

    doc <- list(content = s, meta = meta, annotations = annotations)
    class(doc) <- c("AnnotatedPlainTextDocument",
                    "PlainTextDocument",
                    "TextDocument")

    doc
}

format.AnnotatedPlainTextDocument <-
function(x, ...)
{        
    annotations <- x$annotations
    c(.format_TextDocument(x),
      sprintf("Annotations:  %d, length(s): %s",
              length(annotations),
              paste(sapply(annotations, length), collapse = "/")),
      sprintf("Content:  chars: %d",
              nchar(x$content)))
}
    
## print.AnnotatedPlainTextDocument <-
## function(x, ...)
## {
##     annotations <- x$annotations
##     writeLines(sprintf("<<AnnotatedPlainTextDocument (annotations: %d, length(s): %s)>>",
##                        length(annotations),
##                        paste(sapply(annotations, length),
##                              collapse = "/")))
##     invisible(x)
## }

content.AnnotatedPlainTextDocument <-
function(x)
    x$content

`content<-.AnnotatedPlainTextDocument` <-
function(x, value)
    stop("content modification is not possible for AnnotatedPlainTextDocument objects")

## meta.AnnotatedPlainTextDocument <-
## function(x, tag = NULL, ...)
##     if(is.null(tag)) x$meta else x$meta[[tag]]

## `meta<-.AnnotatedPlainTextDocument` <-
## function(x, tag = NULL, ..., value)
## {
##     if(is.null(tag))
##         x$meta <- value
##     else
##         x$meta[[tag]] <- value
##     x
## }

as.character.AnnotatedPlainTextDocument <-
function(x, ...)
    x$content

annotations <-
function(x)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    x$annotations
}

## NLTK style functions for high level access 

words.AnnotatedPlainTextDocument <-
function(x, which = 1L, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    s <- x$content
    a <- annotations(x)[[which]]
    ## Could check for word token annotations ...
    s[a[a$type == "word"]]
}

sents.AnnotatedPlainTextDocument <-
function(x, which = 1L, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    s <- x$content
    a <- annotations(x)[[which]]
    .sents_from_annotation_and_text(a, s)
}

.sents_from_annotation_and_text <-
function(a, s)
{
    ## Could check for sentence and word token annotations ...
    s[annotations_in_spans(a[a$type == "word"],
                           a[a$type == "sentence"])]
}

paras.AnnotatedPlainTextDocument <-
function(x, which = 1L, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    s <- x$content
    a <- annotations(x)[[which]]
    ## Could check for paragraph annotations ...
    lapply(annotations_in_spans(a, a[a$type == "paragraph"]),
           .sents_from_annotation_and_text, s)
}

tagged_words.AnnotatedPlainTextDocument <-
function(x, which = 1L, map = NULL, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    s <- x$content
    a <- annotations(x)[[which]]
    ## Could check for word token annotations ...
    a <- a[a$type == "word"]
    if(!is.null(map))
        a <- .map_POS_tags_Annotation(a, map)
    .tagged_words_from_annotation_and_text(a, s)
}

.tagged_words_from_annotation_and_text <-
function(a, s)
{
    ## Could check for POS tag features ...
    pos <- sapply(a$features, `[[`, "POS")
    Tagged_Token(s[a], pos)
}

tagged_sents.AnnotatedPlainTextDocument <-
function(x, which = 1L, map = NULL, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    s <- x$content
    a <- annotations(x)[[which]]
    if(!is.null(map))
        a <- .map_POS_tags_Annotation(a, map)
    .tagged_sents_from_annotation_and_text(a, s)
}

.tagged_sents_from_annotation_and_text <-
function(a, s)
{
    ## Could check for word and sentence token annotations ...
    lapply(annotations_in_spans(a[a$type == "word"],
                                a[a$type == "sentence"]),
           .tagged_words_from_annotation_and_text, s)
}

tagged_paras.AnnotatedPlainTextDocument <-
function(x, which = 1L, map = NULL, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    s <- x$content
    a <- annotations(x)[[which]]
    if(!is.null(map))
        a <- .map_POS_tags_Annotation(a, map)
    ## Could check for paragraph annotations ...
    lapply(annotations_in_spans(a, a[a$type == "paragraph"]),
           .tagged_sents_from_annotation_and_text, s)
}
           
parsed_sents.AnnotatedPlainTextDocument <-
function(x, which = 1L, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    a <- annotations(x)[[which]]
    .parsed_sents_from_annotation(a)
}

.parsed_sents_from_annotation <-
function(a)
{
    ## Could check for sentence token annotations and parse features ...
    ptexts <- sapply(a[a$type == "sentence"]$features, `[[`, "parse")
    lapply(ptexts, Tree_parse)
}

parsed_paras.AnnotatedPlainTextDocument <-
function(x, which = 1L, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    a <- annotations(x)[[which]]
    ## Could check for paragraph annotations ...
    lapply(annotations_in_spans(a, a[a$type == "paragraph"]),
           .parsed_sents_from_annotation)
}

chunked_sents.AnnotatedPlainTextDocument <-
function(x, which = 1L, ...)
{
    if(!inherits(x, "AnnotatedPlainTextDocument"))
        stop("argument 'x' must be an AnnotatedPlainTextDocument object")
    s <- x$content
    a <- annotations(x)[[which]]

    ## Require annotations with POS and chunk_tag features, as obtained
    ## e.g. with the Apache OpenNLP POS tag and chunk annotators.  We
    ## could alternatively use annotations with parse features and
    ## flatten the parse trees.

    ## Could check for word and sentence token annotations ...
    lapply(annotations_in_spans(a[a$type == "word"],
                                a[a$type == "sentence"]),
           function(a) {
               ## Could check for POS and chunk tag features ...
               ptags <- sapply(a$features, `[[`, "POS")
               ctags <- sapply(a$features, `[[`, "chunk_tag")
               words <- s[a]
               chunk_tree_from_chunk_info(words, ptags, ctags)
           })
}

.map_POS_tags_Annotation <-
function(x, map)
{
    map <- POS_tag_mapper(map, meta(x, "POS_tagset"))
    x$features <-
        lapply(x$features,
               function(e) {
                   if(!is.null(pos <- e$POS))
                       e$POS <- map(pos)
                   e
               })
    x
}
