calculate.quantiles <-
function (differences, remove.ratio) 
{
    quantiles <- seq(remove.ratio, 1 - remove.ratio, 0.01)
    differences.sorted <- sort(differences)
    indices <- floor(quantiles * length(differences.sorted)) + 
        1
    indices[length(indices)] = indices[length(indices) - 1]
    differences.quantile <- differences.sorted[indices]
    return(differences.quantile)
}
