stat <-
function(x,k)
{
    x <- as.matrix(x)
    n <- nrow(x)
    index <- 1:k
    summ <- 0

    for (i in k:n)
    {
    summ <- summ + sum(x[index])^2
    index <- index+1
    }

    vr1 <- sum(x^2)/n
    vr2 <- summ/(n*k)

    vr <- vr2/vr1
    tem1 <- 2*(2*k-1)*(k-1)
    tem2 <- 3*k*n

    vrstat <- (vr-1)/sqrt( tem1/tem2 )
}
