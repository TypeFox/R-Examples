# example1b.R -- version 2010-12-27
penalty <- function(param, data) {
    y <- data$model(param, data$tm)
    aux <- abs(y - abs(y))
    sum(aux) * data$ww
}
OF <- function(param, data) {		
    y   <- data$model(param, data$tm)
    aux <- y - data$yM
    res <- max(abs(aux))
    # compute the penalty
    aux <- y - abs(y) # ... aux == zero for nonnegative y 
    aux <- -sum(aux) * data$ww
    res + aux
}