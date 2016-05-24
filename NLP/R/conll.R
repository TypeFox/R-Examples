CoNLLTextDocument <-
function(con, encoding = "unknown", meta = list())
{
    ## <NOTE>
    ## Could make the fields controllable, e.g.
    ##   CoNLLDocument(con, fields = c("word", "POS", "chunk_tag"))
    ## and use "" for something to be skipped.
    ## </NOTE>

    records <- scan(con, what = list("", "", ""), encoding = encoding,
                    quote = NULL, quiet = TRUE, fill = TRUE,
                    blank.lines.skip = FALSE)
    words <- records[[1L]]
    ind <- words == ""

    doc <- list(content =
                data.frame(sent = cumsum(ind) + 1L,
                           word = words,
                           POS = records[[2L]],
                           chunk_tag = records[[3L]],
                           stringsAsFactors = FALSE)[!ind, ],
                meta = meta)
    class(doc) <- c("CoNLLTextDocument", "TextDocument")
    doc
}

format.CoNLLTextDocument <-
function(x, ...)
{
    content <- x$content
    nr <- NROW(content)
    c(.format_TextDocument(x),
      sprintf("Content:  words: %d, sents: %d",
              nr,
              content[[nr, "sent"]]))
}

## print.CoNLLTextDocument <-
## function(x, ...)
## {
##     content <- x$content
##     nr <- NROW(content)
##     writeLines(sprintf("<<CoNLLTextDocument (words: %d, sents: %d)>>",
##                        nr, content[[nr, "sent"]]))
##     invisible(x)
## }

content.CoNLLTextDocument <-
function(x)
    x$content

## meta.CoNLLTextDocument <-
## function(x, tag = NULL, ...)
##     if(is.null(tag)) x$meta else x$meta[[tag]]

## `meta<-.CoNLLTextDocument` <-
## function(x, tag = NULL, ..., value)
## {
##     if(is.null(tag))
##         x$meta <- value
##     else
##         x$meta[[tag]] <- value
##     x
## }

as.character.CoNLLTextDocument <-
words.CoNLLTextDocument <-
function(x, ...)
{
    x$content$word
}

sents.CoNLLTextDocument <-
function(x, ...)
{
    split(x$content$word,
          x$content$sent)
}

tagged_words.CoNLLTextDocument <-
function(x, map = NULL, ...)
{
    if(!is.null(map))
        x <- .map_POS_tags_CoNLLTextDocument(x, map)
    Tagged_Token(x$content$word, x$content$POS)
}

tagged_sents.CoNLLTextDocument <-
function(x, map = NULL, ...)
{
    if(!is.null(map))
        x <- .map_POS_tags_CoNLLTextDocument(x, map)
    split(Tagged_Token(x$content$word, x$content$POS),
          x$content$sent)
}

chunked_sents.CoNLLTextDocument <-
function(x, ...)
{
    Map(chunk_tree_from_chunk_info,
        split(x$content$word, x$content$sent),
        split(x$content$POS, x$content$sent),
        split(x$content$chunk_tag, x$content$sent))
}

.map_POS_tags_CoNLLTextDocument <-
function(x, map)
{
    map <- POS_tag_mapper(map, meta(x, "POS_tagset"))
    x$content$POS <- map(x$content$POS)
    x
}
