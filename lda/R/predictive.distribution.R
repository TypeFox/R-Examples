predictive.distribution <-
function (document_sums, topics, alpha, eta) 
{
    smoothed.topics <- (topics + eta)/rowSums(topics + eta)
    apply(document_sums, 2, function(x) {
        props <- (x + alpha)/sum(x + alpha)
        colSums(props * smoothed.topics)
    })
}
