mmsb.collapsed.gibbs.sampler <-
function (network, K, num.iterations, alpha, beta.prior, initial = NULL, 
    burnin = NULL, trace = 0L) 
{
    beta.prior <- list(rep(as.double(beta.prior[[1]]), length.out = K * 
        K), rep(as.double(beta.prior[[2]]), length.out = K * 
        K))
    empty.docs <- rep(list(matrix(integer(0), nrow = 2, ncol = 0)), 
        dim(network)[1])
    retval <- .pairwise.link.lda.collapsed.gibbs.sampler(empty.docs, 
        K, NULL, num.iterations, alpha, 1, beta.prior, network, 
        initial, burnin, trace = trace)
    retval[c("document_sums", "net.assignments.left", "net.assignments.right", 
        "blocks.neg", "blocks.pos", if (is.null(burnin)) NULL else "document_expects")]
}
