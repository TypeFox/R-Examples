## Authors: Ingo Feinerer, Kurt Hornik

TermDocumentMatrix_classes <-
    c("TermDocumentMatrix", "simple_triplet_matrix")
DocumentTermMatrix_classes <-
    c("DocumentTermMatrix", "simple_triplet_matrix")

.TermDocumentMatrix <-
function(x, weighting)
{
    x <- as.simple_triplet_matrix(x)
    if(!is.null(dimnames(x)))
        names(dimnames(x)) <- c("Terms", "Docs")
    class(x) <- TermDocumentMatrix_classes
    ## <NOTE>
    ## Note that if weighting is a weight function, it already needs to
    ## know whether we have a term-document or document-term matrix.
    ##
    ## Ideally we would require weighting to be a WeightFunction object
    ## or a character string of length 2.  But then
    ##   dtm <- DocumentTermMatrix(crude,
    ##                             control = list(weighting =
    ##                                            function(x)
    ##                                            weightTfIdf(x, normalize =
    ##                                                        FALSE),
    ##                                            stopwords = TRUE))
    ## in example("DocumentTermMatrix") fails [because weightTfIdf() is
    ## a weight function and not a weight function generator ...]
    ## Hence, for now, instead of
    ##   if(inherits(weighting, "WeightFunction"))
    ##      x <- weighting(x)
    ## use
    if(is.function(weighting))
        x <- weighting(x)
    ## and hope for the best ...
    ## </NOTE>
    else if(is.character(weighting) && (length(weighting) == 2L))
        attr(x, "weighting") <- weighting
    else
        stop("invalid weighting")
    x
}

TermDocumentMatrix <-
function(x, control = list())
    UseMethod("TermDocumentMatrix", x)

TermDocumentMatrix.PCorpus <-
TermDocumentMatrix.VCorpus <-
function(x, control = list())
{
    stopifnot(is.list(control))

    tflist <- mclapply(unname(content(x)), termFreq, control)
    tflist <- lapply(tflist, function(y) y[y > 0])

    v <- unlist(tflist)
    i <- names(v)
    allTerms <- sort(unique(as.character(if (is.null(control$dictionary)) i
                                         else control$dictionary)))
    i <- match(i, allTerms)
    j <- rep(seq_along(x), sapply(tflist, length))
    docs <- as.character(meta(x, "id", "local"))
    if (length(docs) != length(x)) {
        warning("invalid document identifiers")
        docs <- NULL
    }

    m <- simple_triplet_matrix(i = i, j = j, v = as.numeric(v),
                               nrow = length(allTerms),
                               ncol = length(x),
                               dimnames =
                                 list(Terms = allTerms,
                                      Docs = docs))

    bg <- control$bounds$global
    if (length(bg) == 2L && is.numeric(bg)) {
        rs <- row_sums(m > 0)
        m <- m[(rs >= bg[1]) & (rs <= bg[2]), ]
    }

    weighting <- control$weighting
    if (is.null(weighting))
        weighting <- weightTf

    .TermDocumentMatrix(m, weighting)
}

DocumentTermMatrix <-
function(x, control = list())
    t(TermDocumentMatrix(x, control))

as.TermDocumentMatrix <-
function(x, ...)
    UseMethod("as.TermDocumentMatrix")
as.TermDocumentMatrix.TermDocumentMatrix <-
function(x, ...)
    x
as.TermDocumentMatrix.DocumentTermMatrix <-
function(x, ...)
    t(x)
as.TermDocumentMatrix.term_frequency <-
as.TermDocumentMatrix.textcnt <-
function(x, ...)
{
    m <- simple_triplet_matrix(i = seq_along(x),
                               j = rep(1, length(x)),
                               v = as.numeric(x),
                               nrow = length(x),
                               ncol = 1,
                               dimnames =
                               list(Terms = names(x),
                                    Docs = NA_character_))

    .TermDocumentMatrix(m, weightTf)
}
as.TermDocumentMatrix.default <-
function(x, weighting, ...)
    .TermDocumentMatrix(x, weighting)

as.DocumentTermMatrix <-
function(x, ...)
    UseMethod("as.DocumentTermMatrix")
as.DocumentTermMatrix.DocumentTermMatrix <-
function(x, ...)
    x
as.DocumentTermMatrix.TermDocumentMatrix <-
function(x, ...)
    t(x)
as.DocumentTermMatrix.term_frequency <-
as.DocumentTermMatrix.textcnt <-
function(x, ...)
    t(as.TermDocumentMatrix(x))
as.DocumentTermMatrix.default <-
function(x, weighting, ...)
{
    x <- as.simple_triplet_matrix(x)
    t(.TermDocumentMatrix(t(x), weighting))
}

t.TermDocumentMatrix <-
t.DocumentTermMatrix <-
function(x)
{
    m <- NextMethod("t")
    attr(m, "weighting") <- attr(x, "weighting")
    class(m) <- if(inherits(x, "DocumentTermMatrix"))
        TermDocumentMatrix_classes
    else
        DocumentTermMatrix_classes
    m
}

termFreq <-
function(doc, control = list())
{
    stopifnot(inherits(doc, "TextDocument"), is.list(control))

    ## Tokenize the corpus
    .tokenize <- control$tokenize
    if (is.null(.tokenize) || identical(.tokenize, "words"))
        .tokenize <- words
    else if (identical(.tokenize, "MC"))
        .tokenize <- MC_tokenizer
    else if (identical(.tokenize, "scan"))
        .tokenize <- scan_tokenizer
    else if (NLP::is.Span_Tokenizer(.tokenize))
        .tokenize <- NLP::as.Token_Tokenizer(.tokenize)
    if (is.function(.tokenize))
        txt <- .tokenize(doc)
    else
        stop("invalid tokenizer")

    ## Conversion to lower characters
    .tolower <- control$tolower
    if (is.null(.tolower) || isTRUE(.tolower))
        .tolower <- tolower
    if (is.function(.tolower))
        txt <- .tolower(txt)

    ## Punctuation removal
    .removePunctuation <- control$removePunctuation
    if (isTRUE(.removePunctuation))
        .removePunctuation <- removePunctuation
    else if (is.list(.removePunctuation))
        .removePunctuation <-
            function(x) do.call(removePunctuation,
                                c(list(x), control$removePunctuation))

    ## Number removal
    .removeNumbers <- control$removeNumbers
    if (isTRUE(.removeNumbers))
        .removeNumbers <- removeNumbers

    ## Stopword filtering
    .stopwords <- control$stopwords
    if (isTRUE(.stopwords))
        .stopwords <- function(x) x[is.na(match(x, stopwords(meta(doc, "language"))))] 
    else if (is.character(.stopwords))
        .stopwords <- function(x) x[is.na(match(x, control$stopwords))]

    ## Stemming
    .stemming <- control$stemming
    if (isTRUE(.stemming))
        .stemming <- function(x) stemDocument(x, meta(doc, "language"))

    ## Default order for options which support reordering
    or <- c("removePunctuation", "removeNumbers", "stopwords", "stemming")

    ## Process control options in specified order
    nc <- names(control)
    n <- nc[nc %in% or]
    for (name in sprintf(".%s", c(n, setdiff(or, n)))) {
        g <- get(name)
        if (is.function(g))
            txt <- g(txt)
    }

    ## Check if the document content is NULL
    if (is.null(txt))
        return(setNames(integer(0), character(0)))

    ## If dictionary is set tabulate against it
    dictionary <- control$dictionary
    tab <-  if (is.null(dictionary))
        table(txt)
    else
        table(factor(txt, levels = dictionary))

    ## Ensure local bounds
    bl <- control$bounds$local
    if (length(bl) == 2L && is.numeric(bl))
        tab <- tab[(tab >= bl[1]) & (tab <= bl[2])]

    ## Filter out too short or too long terms
    nc <- nchar(names(tab), type = "chars")
    wl <- control$wordLengths
    lb <- if (is.numeric(wl[1])) wl[1] else 3
    ub <- if (is.numeric(wl[2])) wl[2] else Inf
    tab <- tab[(nc >= lb) & (nc <= ub)]

    ## Return named integer
    storage.mode(tab) <- "integer"
    class(tab) <- c("term_frequency", class(tab))
    tab
}

print.TermDocumentMatrix <-
print.DocumentTermMatrix <-
function(x, ...)
{
    format <- c("term", "document")
    if (inherits(x, "DocumentTermMatrix"))
        format <- rev(format)
    writeLines(sprintf("<<%s (%ss: %d, %ss: %d)>>",
                       class(x)[1], format[1L], nrow(x), format[2L], ncol(x)))
    writeLines(sprintf("Non-/sparse entries: %d/%.0f",
                length(x$v), prod(dim(x)) - length(x$v)))
    sparsity <- if (!prod(dim(x))) 100
        else round((1 - length(x$v)/prod(dim(x))) * 100)
    writeLines(sprintf("Sparsity           : %s%%", sparsity))
    writeLines(sprintf("Maximal term length: %s",
                       max(nchar(Terms(x), type = "chars"), 0)))
    writeLines(sprintf("Weighting          : %s (%s)",
                       attr(x, "weighting")[1L], attr(x, "weighting")[2L]))
    invisible(x)
}

inspect.TermDocumentMatrix <-
inspect.DocumentTermMatrix <-
function(x)
{
    print(x)
    cat("\n")
    print(as.matrix(x))
}

`[.TermDocumentMatrix` <-
`[.DocumentTermMatrix` <-
function(x, i, j, ..., drop)
{
    m <- NextMethod("[")
    attr(m, "weighting") <- attr(x, "weighting")
    class(m) <- if (inherits(x, "DocumentTermMatrix"))
        DocumentTermMatrix_classes
    else
        TermDocumentMatrix_classes
    m
}

`dimnames<-.DocumentTermMatrix` <-
function(x, value)
{
    x <- NextMethod("dimnames<-")
    dnx <- x$dimnames
    if(!is.null(dnx))
        names(dnx) <- c("Docs", "Terms")
    x$dimnames <- dnx
    x
}

`dimnames<-.TermDocumentMatrix` <-
function(x, value)
{
    x <- NextMethod("dimnames<-")
    dnx <- x$dimnames
    if(!is.null(dnx))
        names(dnx) <- c("Terms", "Docs")
    x$dimnames <- dnx
    x
}

nDocs <-
function(x)
    UseMethod("nDocs")

nTerms <-
function(x)
    UseMethod("nTerms")

nDocs.DocumentTermMatrix <-
nTerms.TermDocumentMatrix <-
function(x)
    x$nrow

nDocs.TermDocumentMatrix <-
nTerms.DocumentTermMatrix <-
function(x)
    x$ncol

Docs <-
function(x)
    UseMethod("Docs")

Terms <-
function(x)
    UseMethod("Terms")

Docs.DocumentTermMatrix <-
Terms.TermDocumentMatrix <-
function(x)
{
    s <- x$dimnames[[1L]]
    if(is.null(s))
        s <- rep.int(NA_character_, x$nrow)
    s
}

Docs.TermDocumentMatrix <-
Terms.DocumentTermMatrix <-
function(x)
{
    s <- x$dimnames[[2L]]
    if(is.null(s))
        s <- rep.int(NA_character_, x$ncol)
    s
}

c.term_frequency <-
function(..., recursive = FALSE)
{
    do.call("c", lapply(list(...), as.TermDocumentMatrix))
}

c.TermDocumentMatrix <-
function(..., recursive = FALSE)
{
    m <- lapply(list(...), as.TermDocumentMatrix)

    if(length(m) == 1L)
        return(m[[1L]])

    weighting <- attr(m[[1L]], "weighting")

    allTermsNonUnique <- unlist(lapply(m, function(x) Terms(x)[x$i]))
    allTerms <- unique(allTermsNonUnique)
    allDocs <- unlist(lapply(m, Docs))

    cs <- cumsum(lapply(m, nDocs))
    cs <- c(0, cs[-length(cs)])
    j <- lapply(m, "[[", "j")

    m <- simple_triplet_matrix(i = match(allTermsNonUnique, allTerms),
                               j = unlist(j) + rep.int(cs, sapply(j, length)),
                               v = unlist(lapply(m, "[[", "v")),
                               nrow = length(allTerms),
                               ncol = length(allDocs),
                               dimnames =
                               list(Terms = allTerms,
                                    Docs = allDocs))
    ## <NOTE>
    ## - We assume that all arguments have the same weighting
    ## - Even if all matrices have the same input weighting it might be necessary
    ##   to take additional steps (e.g., normalization for tf-idf or check for
    ##   (0,1)-range for binary tf)
    ## </NOTE>
    .TermDocumentMatrix(m, weighting)
}

c.DocumentTermMatrix <-
function(..., recursive = FALSE)
{
    t(do.call("c", lapply(list(...), as.TermDocumentMatrix)))
}

findFreqTerms <-
function(x, lowfreq = 0, highfreq = Inf)
{
    stopifnot(inherits(x, c("DocumentTermMatrix", "TermDocumentMatrix")),
              is.numeric(lowfreq), is.numeric(highfreq))

    if (inherits(x, "DocumentTermMatrix")) x <- t(x)
    rs <- slam::row_sums(x)
    names(rs[rs >= lowfreq & rs <= highfreq])
}

findAssocs <-
function(x, terms, corlimit)
    UseMethod("findAssocs", x)
findAssocs.TermDocumentMatrix <-
function(x, terms, corlimit)
    findAssocs(t(x), terms, corlimit)
findAssocs.DocumentTermMatrix <-
function(x, terms, corlimit)
{
    stopifnot(is.character(terms), is.numeric(corlimit),
              corlimit >= 0, corlimit <= 1)

    j <- match(unique(terms), Terms(x), nomatch = 0L)
    suppressWarnings(
        findAssocs(slam::crossapply_simple_triplet_matrix(x[, j], x[, -j], cor),
                   terms, rep_len(corlimit, length(terms))))
}
findAssocs.matrix <-
function(x, terms, corlimit)
{
    stopifnot(is.numeric(x))

    i <- match(terms, rownames(x), nomatch = 0L)
    names(i) <- terms
    Map(function(i, cl) {
        xi <- x[i, ]
        t <- sort(round(xi[which(xi >= cl)], 2), TRUE)
        if (!length(t))
            names(t) <- NULL
        t
        }, i, corlimit)
}

removeSparseTerms <-
function(x, sparse)
{
    stopifnot(inherits(x, c("DocumentTermMatrix", "TermDocumentMatrix")),
              is.numeric(sparse), sparse > 0, sparse < 1)

    m <- if (inherits(x, "DocumentTermMatrix")) t(x) else x
    t <- table(m$i) > m$ncol * (1 - sparse)
    termIndex <- as.numeric(names(t[t]))
    if (inherits(x, "DocumentTermMatrix")) x[, termIndex] else x[termIndex,]
}

CategorizedDocumentTermMatrix <-
function(x, c)
{
    if(inherits(x, "TermDocumentMatrix"))
        x <- t(x)
    else if(!inherits(x, "DocumentTermMatrix"))
        stop("wrong class")

    if(length(c) != nDocs(x))
        stop("invalid category ids")

    attr(x, "Category") <- c

    class(x) <- c("CategorizedDocumentTermMatrix",
                  DocumentTermMatrix_classes)

    x
}
