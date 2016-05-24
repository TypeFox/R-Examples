# Predict the category of a new observation
# @param nn A probabilistic neural network already trained
# @param X A vector describing a new observation
guess.category <- function(nn, X) {
    probs <- guess.probabilities.of.each.category(nn, X)
    if(is.na(probs[1])) return(NA)
    category <- names(probs[probs == max(probs)])
    return(category)
}
