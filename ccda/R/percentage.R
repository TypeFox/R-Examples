percentage <-
function (dataset, starting_vector, prior) 
{
    if (prior == "proportions") {
        confmtx = table(predict(lda(dataset, starting_vector), 
            dataset)$class, true = starting_vector)
    }
    if (prior == "equal") {
        confmtx = table(predict(lda(dataset, starting_vector, prior = c(rep(1/length(table(starting_vector)), 
            times = length(table(starting_vector))))), dataset)$class, 
            true = starting_vector)
    }
    perctg = sum(diag(confmtx))/sum(confmtx)
    return(perctg)
}
