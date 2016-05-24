`rgg` <-
function(m, alpha, beta, n)
{
# Generation of a m Gamma-gamma variables i.i.d

        l <- 1/beta * rgamma(m, alpha)
        1/l * rgamma(m, n)
}

