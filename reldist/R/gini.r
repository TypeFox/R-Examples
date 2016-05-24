gini <- function(x, weights=rep(1,length=length(x))){
        ox <- order(x)
        x <- x[ox]
        weights <- weights[ox]/sum(weights)
        p <- cumsum(weights)
        nu <- cumsum(weights*x)
        n <- length(nu)
        nu <- nu / nu[n]
        sum(nu[-1]*p[-n]) - sum(nu[-n]*p[-1])
}
