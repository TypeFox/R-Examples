top.topic.documents <-
function (document_sums, num.documents = 20, alpha = 0.1) 
{
    normalized.sums <- t(document_sums + alpha)/colSums(document_sums + 
        alpha)
    apply(normalized.sums, 2, function(x) order(x, decreasing = TRUE)[1:num.documents])
}
