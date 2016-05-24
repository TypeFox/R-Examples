mode.lands.t <-
function (nu, zeta) 
{
    if (nu < 2) 
        stop("'nu' must be greater than or equal to 2")
    fcn <- function(theta, nu, zeta) {
        -(nu - 1) * tan(theta) + (1 + nu) * zeta * cos(theta)
    }
    theta.max <- uniroot(fcn, interval = c(-pi/2, pi/2), nu = nu, 
        zeta = zeta)$root
    sqrt(nu) * tan(theta.max)
}
