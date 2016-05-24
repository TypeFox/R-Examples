inner.prod.iter <-
function(x) {
n <- length(x)
I.plus <- I.minus <- I.prod <- rep(0, n-1)
I.plus[1] <- sqrt(1 - 1/n) * x[1]
I.minus[1] <- 1/sqrt(n^2 - n) * sum(x[2:n])

if (n-2) for (m in 1:(n-2)) {
factor <- sqrt((n-m-1)*m/(m+1)/(n-m))
I.plus[m+1] <- I.plus[m] * factor + x[m+1] * sqrt(1/(m+1) - 1/n)
I.minus[m+1] <- I.minus[m] / factor - x[m+1] / sqrt(n^2/(m+1)-n)
}

I.prod <- I.plus - I.minus

return(I.prod)
}

