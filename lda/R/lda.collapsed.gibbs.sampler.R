lda.collapsed.gibbs.sampler <-
function (documents, K, vocab, num.iterations, alpha, eta, initial = NULL, 
    burnin = NULL, compute.log.likelihood = FALSE, trace = 0L, 
    freeze.topics = FALSE) 
{
    if (class(vocab) == "list") {
        lengths <- as.integer(sapply(vocab, length))
        all.vocab <- do.call(c, vocab)
    }
    else {
        lengths <- as.integer(length(vocab))
        all.vocab <- vocab
    }
    retval <- structure(.Call("collapsedGibbsSampler", documents, 
        as.integer(K), lengths, as.integer(num.iterations), as.double(alpha), 
        as.double(eta), NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
        initial, as.integer(burnin), as.logical(compute.log.likelihood), 
        trace, as.logical(freeze.topics)), names = c("assignments", 
        "topics", "topic_sums", "document_sums", if (is.null(burnin)) NA else "document_expects", 
        NA, NA, NA, NA, if (compute.log.likelihood) "log.likelihoods" else NA))
    colnames(retval$topics) <- all.vocab
    retval
}
