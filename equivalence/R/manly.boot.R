"manly.boot" <-
function (x, i, Epsilon=1)
{
    reps <- 50
    x <- x[i]
    mean <- mean(x)
    n <- length(x)

    pmanly1.boot <- function(d, i)
    {  randomize <- 1 - 2*rbinom(n, 1, 0.5)
       d1.boot <- (d[i]+Epsilon)*randomize
       return(mean(d1.boot) >= mean+Epsilon) }
    pmanly2.boot <- function(d, i)
    {  randomize <- 1 - 2*rbinom(n, 1, 0.5)
       d2.boot <- (d[i]-Epsilon)*randomize
       return(mean(d2.boot) <= mean-Epsilon) }

    pvalue <- mean(boot(x, pmanly1.boot, reps)$t) +
              mean(boot(x, pmanly2.boot, reps)$t)

    result <- (pvalue < 0.05)
}

