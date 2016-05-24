rtm.collapsed.gibbs.sampler <-
function (documents, links, K, vocab, num.iterations, alpha, 
    eta, beta, trace = 0L, test.start = length(documents) + 1L) 
{
    retval <- structure(.Call("rtm", documents, links, as.integer(K), 
        length(vocab), as.integer(num.iterations), as.double(alpha), 
        as.double(eta), rep(as.double(beta), length.out = K), 
        trace, test.start), names = c("assignments", "topics", 
        "topic_sums", "document_sums", "log.likelihoods"))
    colnames(retval$topics) <- vocab
    retval
}
