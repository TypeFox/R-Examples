calculate.absmax <-
function (differences.quantile) 
{
    differences.quantile.abs <- abs(differences.quantile)
    differences.quantile.abs[differences.quantile.abs == Inf] <- 0
    absmax <- max(differences.quantile.abs)
    return(absmax)
}
