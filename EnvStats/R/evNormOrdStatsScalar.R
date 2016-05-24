evNormOrdStatsScalar <-
function (r = 1, n = 1, approximate = FALSE) 
{
    if (any(length.list(r, n) != 1) || r != trunc(r) || n != 
        trunc(n) || r < 1 || n < 1 || r > n) 
        stop(paste("'r' and 'n' must be positive integers,", 
            "and 'r' must be between 1 and 'n'"))
    if (approximate) 
        ev <- qnorm((r - 0.375)/(n + 0.25))
    else {
        if (is.odd(n) && r == (n + 1)/2) 
            ev <- 0
        else {
            integrand.fcn <- function(x, r.weird, n.weird) {
                x * (pnorm(x)^(r.weird - 1)) * ((1 - pnorm(x))^(n.weird - 
                  r.weird)) * dnorm(x)
            }
            ev <- exp(log(r) + lchoose(n, r)) * integrate(integrand.fcn, 
                lower = -Inf, upper = Inf, r.weird = r, n.weird = n)$value
        }
    }
    ev
}
