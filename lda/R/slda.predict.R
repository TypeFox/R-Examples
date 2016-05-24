slda.predict <-
function (documents, topics, model, alpha, eta, num.iterations = 100, 
    average.iterations = 50, trace = 0L) 
{
    doc_sums_count <- slda.predict.docsums(documents, topics, 
        alpha, eta, num.iterations, average.iterations, trace)
    props <- t(doc_sums_count)/colSums(doc_sums_count)
     props %*% matrix(coef(model),nrow=ncol(props), byrow=T)
}
