R2.global <-
function(z, partition)
{
    n <- tapply(partition, partition, length)
    m <- tapply(z, partition, mean)

    mT <- mean(z)
    v <- var(z)

    if (var(z) == 0)
    {
        return(0)
    }

    R2 <- sum(n*(m-mT)**2) / (v*(length(z)-1))

    return(R2)
}
