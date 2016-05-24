lda.cvb0 <- function(documents, K, vocab, num.iterations, alpha, eta, trace = 0L) {
    retval <- structure(
      .Call(
        "cvb0", 
        documents, 
        as.integer(K), 
        as.integer(length(vocab)), 
        as.integer(num.iterations), 
        as.double(alpha), 
        as.double(eta),
        as.integer(trace)
      ), 
      names = c(
        "assignments", 
        "topics", 
        "topic_sums", 
        "document_sums"
      )
    )
    colnames(retval$topics) <- vocab
    retval
}
