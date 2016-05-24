getmeans <-
function(data, factors)
{
    groups <- levels(as.factor(factors))
    means <- array(NA, c(length(groups), dim(data)[2]))
    for(i in 1:length(groups))
        means[i,] <- apply(data[factors == groups[i], , drop = FALSE], 2, mean)
    colnames(means) <- colnames(data)
    return(means)
}
