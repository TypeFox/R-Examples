predictive.link.probability <-
function (edgelist, document_sums, alpha, beta) 
{
    beta <- rep(beta, length.out = dim(document_sums)[1])
    apply(edgelist + 1, 1, function(x) {
        exp((beta * (document_sums[, x[1]] + alpha)) %*% (document_sums[, 
            x[2]] + alpha)/sum(document_sums[, x[1]] + alpha)/sum(document_sums[, 
            x[2]] + alpha) - max(beta, 0))
    })
}
