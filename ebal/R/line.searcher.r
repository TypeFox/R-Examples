# function to conduct line search for optimal step length
line.searcher <- function(
                        Base.weight,
                        Co.x,
                        Tr.total,
                        coefs,
                        Newton,
                        ss)
 {
    weights.temp <- c(exp(Co.x %*% (coefs - (ss * Newton) )))
    #weights.temp[is.infinite(weights.temp)] <- 100
    weights.temp <- weights.temp * Base.weight
    Co.x.agg     <- c(weights.temp %*% Co.x)
    maxdiff      <- max(abs(Co.x.agg-Tr.total))
    return(maxdiff)
}