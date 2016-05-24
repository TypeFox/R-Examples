## Tokenizers break text up into words, phrases, symbols, or other
## meaningful elements called tokens, see e.g.
## <http://en.wikipedia.org/wiki/Tokenization_%28lexical_analysis%29>.
## This can be accomplished by returning the sequence of tokens, or the
## corresponding spans (character start and end positions).
## Apache OpenNLP provides a Tokenizer interface, with methods
## String[] tokenize() and Span[] tokenizePos() for the two variants.
## See e.g.
## <https://opennlp.apache.org/documentation/1.5.3/apidocs/opennlp-tools/opennlp/tools/tokenize/Tokenizer.html>.
## NLTK provides an interface class nltk.tokenize.api.TokenizerI, for
## which subclasses must define a tokenize() method, and can define a
## span_tokenize() method.
## See e.g. <http://www.nltk.org/api/nltk.tokenize.html>.
## In R, this could be mimicked by having two generics for getting the
## tokens or spans, and have a virtual Tokenizer class for which
## extension classes must provide methods for at least one of the
## generics.
## However, it seems more natural to have tokenizers be *functions*
## (instead of interface classes) which can be called directly (instead
## of calling the respective generics), and have two "kinds" of such
## functions: token tokenizers and span tokenizers.  We use the class
## information to indicate the kind, which in turn allows to provide a
## generic mechanism for mapping between the two kinds (straightforward
## when going from spans to tokens, doable for the opposite direction).
## This also allows to "extract" both kinds of tokenizers from suitable
## annotators or annotator pipelines.
## For now, there is no underlying virtual Tokenizer class.

### * Span tokenizers

Span_Tokenizer <-
function(f, meta = list())
{
    attr(f, "meta") <- meta
    class(f) <- "Span_Tokenizer"
    f
}

as.Span_Tokenizer <-
function(x, ...)
    UseMethod("as.Span_Tokenizer")

as.Span_Tokenizer.Span_Tokenizer <-
function(x, ...)
    x

## For now, pass metadata as is.
as.Span_Tokenizer.Token_Tokenizer <-
function(x, ...)
{
    f <- function(s) {
        s <- as.String(s)
        spans_from_tokens(s, x(s))
    }
    Span_Tokenizer(f, meta(x))
}

## For now, do not pass metadata.
as.Span_Tokenizer.Annotator <-
as.Span_Tokenizer.Annotator_Pipeline <-
function(x, type = "word", ...)
{
    f <- function(s) {
        a <- x(as.String(s))
        as.Span(a[a$type == "word", ])
    }
    Span_Tokenizer(f)
}
        
is.Span_Tokenizer <-
function(x)
    inherits(x, "Span_Tokenizer")

format.Span_Tokenizer <-
function(x, ...)
{
    d <- meta(x, "description")
    if(is.null(d)) {
        "A span tokenizer."
    } else {
        c("A span tokenizer, with description",
          strwrap(d, indent = 2L, exdent = 2L))
    }
}

### * Token tokenizers

Token_Tokenizer <-
function(f, meta = list())    
{
    attr(f, "meta") <- meta
    class(f) <- "Token_Tokenizer"
    f
}

as.Token_Tokenizer <-
function(x, ...)
    UseMethod("as.Token_Tokenizer")

as.Token_Tokenizer.Token_Tokenizer <-
function(x, ...)
    x

## For now, pass metadata as is.
as.Token_Tokenizer.Span_Tokenizer <-
function(x, ...)
{
    f <- function(s) {
        s <- as.String(s)
        s[x(s)]
    }
    Token_Tokenizer(f, meta(x))
}
    
## For now, do not pass metadata.
as.Token_Tokenizer.Annotator <-
as.Token_Tokenizer.Annotator_Pipeline <-
function(x, type = "word", ...)
{
    f <- function(s) {
        s <- as.String(s)
        a <- x(s)
        s[a[a$type == "word", ]]
    }
    Token_Tokenizer(f)
}

is.Token_Tokenizer <-
function(x)
    inherits(x, "Token_Tokenizer")

format.Token_Tokenizer <-
function(x, ...)
{
    d <- meta(x, "description")
    if(is.null(d)) {
        "A token tokenizer."
    } else {
        c("A token tokenizer, with description",
          strwrap(d, indent = 2L, exdent = 2L))
    }
}   

### Regexp span tokenizers a la NLTK.

Regexp_Tokenizer <-
function(pattern, invert = FALSE, ..., meta = list())
{
    force(pattern)
    args <- list(...)
    
    f <- if(invert) {
        ## Pattern gives the separators.
        function(s) {
            s <- as.String(s)
            if(is.na(s) || !nchar(s))
                stop("Need a non-empty string.")
            m <- do.call(gregexpr,
                         c(list(pattern = pattern, text = s), args))[[1L]]
            if((length(m) == 1L) && (m == -1L))
                return(Span(1L, nchar(s)))
            start <- c(1L, m + attr(m, "match.length"))
            end <- c(m - 1L, nchar(s))
            ind <- start <= end
            Span(start[ind], end[ind])
        }
    } else {
        ## Pattern gives the tokens.
        function(s) {
            s <- as.String(s)
            if(is.na(s) || !nchar(s))
                stop("Need a non-empty string.")
            m <- do.call(gregexpr,
                         c(list(pattern = pattern, text = s), args))[[1L]]
            Span(m, m + attr(m, "match.length") - 1L)
        }
    }

    Span_Tokenizer(f, meta)
}

whitespace_tokenizer <-
    Regexp_Tokenizer("\\s+",
                     invert = TRUE,
                     meta = list(description = "Divides strings into substrings by treating any sequence of whitespace characters as a separator."))

blankline_tokenizer <-
    Regexp_Tokenizer("\\s*\n\\s*\\n\\s*",
                     invert = TRUE,
                     meta = list(description = "Divides strings into substrings by treating any sequence of blank lines as a separator."))

wordpunct_tokenizer <-
    Regexp_Tokenizer("\\w+|[^\\w\\s]+",
                     perl = TRUE,
                     meta = list(description = "Divides strings into substrings of alphabetic and (non-whitespace) non-alphabetic characters."))

### * Utilities

spans_from_tokens <-
function(x, tokens)
{
    start <- end <- integer(length(tokens))
    off <- 0L
    for(i in seq_along(tokens)) {
        m <- regexpr(tokens[i], x, fixed = TRUE)
        pos <- m + attr(m, "match.length")
        x <- substring(x, pos)
        start[i] <- off + m
        end[i] <- off <- off + pos - 1L
        
    }
    Span(start, end)
}
