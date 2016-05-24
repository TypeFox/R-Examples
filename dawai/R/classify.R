classify <-
function(y, means, variances, prior, numgroups, dimension)
{
    values <- array(NA, c(dim(y)[1], numgroups))
    if(dimension == 1)
        for(i in 1:numgroups)
            values[, i] <- dnorm(y, mean = means[i], sd = sqrt(variances[, , i]), log = TRUE)
    else
        for(i in 1:numgroups)
            values[, i] <- c(apply(y, 1, dmvnorm, mean = means[{dimension * {i - 1} + 1}:{dimension * i}],
                                   sigma = variances[, , i], log = TRUE))
    values <- t(t(exp(values - apply(values, 1, max)))*c(prior))
    maxClass <- {values == apply(values, 1, max)} * rep(1:numgroups, each = length(as.matrix(y)) / dimension)
    group <- apply(maxClass, 1, function(x) x[x != 0][sample(length(x[x != 0]), 1)])
    return(group)
}
