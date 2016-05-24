getvariances <-
function(data, factors)
{
    dimension <- dim(data)[2]
    vars <- do.call(rbind, lapply(split(as.data.frame(data), as.data.frame(factors)), var))
    variances <- array(t(vars), c(dimension, dimension, length(levels(as.factor(factors)))))
    return(variances)
}