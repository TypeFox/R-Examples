# Predict the probabilities of each category given X
# @param nn An already trained probabilistic neural network
# @param X Pattern from which we have to decide a category. It is a set of measurements represented by a p-dimensional vector
guess.probabilities.of.each.category <- function(nn, X) {
    results <- vector()
    for(category in nn$categories) {
        Xa <- nn$set[nn$set[,nn$category.column] == category,]
        Xa <- as.matrix(Xa[,-nn$category.column])
        results <- c(results, fA(Xa, X, nn$sigma))
    }
    probs <- results / sum(results)
    names(probs) <- nn$categories
    return(probs)
}
