nubbi.collapsed.gibbs.sampler <-
function (contexts, pair.contexts, pairs, K.individual, K.pair, 
    vocab, num.iterations, alpha, eta, xi) 
{
    pair.sources <- cbind(pairs, 0:(length(pair.contexts) - 1) + 
        length(contexts))
    sources <- c(as.list(0:(length(contexts) - 1)), mapply(c, 
        pair.sources[, 1], pair.sources[, 2], pair.sources[, 
            3], SIMPLIFY = FALSE))
    retval <- structure(.Call("nubbi", c(contexts, pair.contexts), 
        sources, rep(0:1, c(length(contexts), length(pair.contexts))), 
        as.integer(c(K.individual, K.pair)), as.integer(length(vocab)), 
        as.integer(num.iterations), as.double(alpha), as.double(eta), 
        as.double(xi)), names = c("topic_assignments", "source_assignments", 
        "topics", "topic_sums", "document_sums", "document_source_sums"))
    retval$source_assignments <- retval$source_assignments[(length(contexts) + 
        1):(length(contexts) + length(pair.contexts))]
    retval$document_source_sums <- matrix(unlist(retval$document_source_sums[(length(contexts) + 
        1):(length(contexts) + length(pair.contexts))]), ncol = 3, 
        byrow = TRUE)
    colnames(retval$topics) <- vocab
    retval
}
