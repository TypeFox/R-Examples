"tost.boot" <-
function (x, i, Epsilon=0.5)
{
    x <- x[i]
    mean <- mean(x)
    std <- sd(x)
    n = length(x)

    tint <- (std/sqrt(n)) * qt(0.95, n-1)
    result <- ((mean + tint) < Epsilon & (mean - tint) > -Epsilon)
}

