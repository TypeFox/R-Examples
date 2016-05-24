lexicalize <-
function (doclines, sep = " ", lower = TRUE, count = 1L, vocab = NULL) 
{
    split.docs <- strsplit(if (lower) 
        tolower(doclines)
    else doclines, sep, fixed = TRUE)
    if (is.null(vocab)) {
        local.vocab <- unique(unlist(split.docs))
    }
    else {
        local.vocab <- vocab
    }
    split.docs <- lapply(split.docs, function(x) {
        m <- match(x, local.vocab)
        m <- m[!is.na(m)]
        rbind(m - 1L, rep(as.integer(count), length(m)))
    })
    if (is.null(vocab)) {
        list(documents = split.docs, vocab = local.vocab)
    }
    else {
        split.docs
    }
}
