ConditionalDistribution <-
function(channel, index = NA, x = NA) {
    ## Returns the conditional distribution P(y | x) for a given channel
    ## and input index.

    if(is.na(index)) {
        index <- which(channel$x == x)[1]
    }
    if(is.na(index)) {
        stop("Invalid index for conditional distribution.")
    }
    
    p <- channel$Aj[[index]] * channel$q
    p <- p / sum(p)
    
    list(y = channel$y, p = p)
}
