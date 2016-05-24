



# add perturbations to weights, without changing the sort order of nonunique 
# elements
#
#
perturbWeights <- function(x){
    d <- diff(sort(x))
    p <- min(x, d[d > 0]) / 100
    x + runif(length(x), max=p)
}

