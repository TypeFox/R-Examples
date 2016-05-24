Dunnetts.K.F1 <-
function (x, k, rho) 
{
    fcn.to.integrate <- function(y, x, k, rho) {
        pnorm((x + sqrt(rho) * y)/sqrt(1 - rho))^k * dnorm(y)
    }
    integrate(fcn.to.integrate, lower = -Inf, upper = Inf, x = x, 
        k = k, rho = rho)$value
}
