# Author: Ingo Feinerer
# Transformations

tm_map <-
function(x, FUN, ...)
    UseMethod("tm_map", x)
tm_map.VCorpus <-
function(x, FUN, ..., lazy = FALSE)
{
    # Lazy mapping
    if (lazy) {
        fun <- function(x) FUN(x, ...)
        if (is.null(x$lazy))
            x$lazy <- list(index = rep(TRUE, length(x)), maps = list(fun))
        else
            x$lazy$maps <- c(x$lazy$maps, list(fun))
    } else
        x$content <- mclapply(content(x), FUN, ...)
    x
}
tm_map.PCorpus <-
function(x, FUN, ...)
{
    db <- filehash::dbInit(x$dbcontrol[["dbName"]], x$dbcontrol[["dbType"]])
    for (i in seq_along(x))
        db[[x$content[[i]]]] <- FUN(x[[i]], ...)
    filehash::dbReorganize(db)
    x
}

# Materialize lazy mappings
materialize <-
function(x, range = seq_along(x))
{
    if (!is.null(x$lazy)) {
       i <- (seq_along(x) %in% range) & x$lazy$index
       if (any(i)) {
           x$content[i] <-
               mclapply(x$content[i], function(d) tm_reduce(d, x$lazy$maps))
           x$lazy$index[i] <- FALSE
       }

       # Clean up if everything is materialized
       if (!any(x$lazy$index))
           x["lazy"] <- list(NULL)
    }
    x
}

tm_reduce <-
function(x, tmFuns, ...)
    Reduce(function(f, ...) f(...), tmFuns, x, right = TRUE)

getTransformations <-
function()
    c("removeNumbers", "removePunctuation", "removeWords", "stemDocument",
      "stripWhitespace")

content_transformer <-
function(FUN)
    function(x, ...) {
        content(x) <- FUN(content(x), ...)
        x
    }

removeNumbers <-
function(x)
    UseMethod("removeNumbers", x)
removeNumbers.character <-
function(x)
    gsub("[[:digit:]]+", "", x)
removeNumbers.PlainTextDocument <-
    content_transformer(removeNumbers.character)

removePunctuation <-
function(x, preserve_intra_word_dashes = FALSE)
    UseMethod("removePunctuation", x)
removePunctuation.character <-
function(x, preserve_intra_word_dashes = FALSE)
{
    if (!preserve_intra_word_dashes)
        gsub("[[:punct:]]+", "", x)
    else {
        # Assume there are no ASCII 1 characters.
        x <- gsub("(\\w)-(\\w)", "\\1\1\\2", x)
        x <- gsub("[[:punct:]]+", "", x)
        gsub("\1", "-", x, fixed = TRUE)
    }
}
removePunctuation.PlainTextDocument <-
    content_transformer(removePunctuation.character)

removeWords <-
function(x, words)
    UseMethod("removeWords", x)
# Improvements by Kurt Hornik
removeWords.character <-
function(x, words)
    gsub(sprintf("(*UCP)\\b(%s)\\b",
                 paste(sort(words, decreasing = TRUE), collapse = "|")),
         "", x, perl = TRUE)
removeWords.PlainTextDocument <-
    content_transformer(removeWords.character)

stemDocument <-
function(x, language = "english")
    UseMethod("stemDocument", x)
stemDocument.character <-
function(x, language = "english")
    SnowballC::wordStem(x, as.character(language))
stemDocument.PlainTextDocument <-
function(x, language = meta(x, "language"))
{
    language <- as.character(language)
    if (identical(language, "") ||
        identical(language, character(0)) ||
        is.na(language))
        language <- "english"

    s <- unlist(lapply(content(x),
      function(x) paste(stemDocument.character(unlist(strsplit(x, "[[:blank:]]")),
                                               language),
                        collapse = " ")))
    content(x) <- if (is.character(s)) s else ""
    x
}

stripWhitespace <-
function(x)
    UseMethod("stripWhitespace", x)
stripWhitespace.character <-
function(x)
    gsub("[[:space:]]+", " ", x)
stripWhitespace.PlainTextDocument <-
    content_transformer(stripWhitespace.character)
