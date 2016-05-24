dlands.t <-
function (tau, nu, zeta) 
{
    if (nu < 2) 
        stop("'nu' must be greater than or equal to 2")
    if (abs(tau) > .Machine$double.xmax) 
        return(0)
    else {
        fcn <- function(theta, nu, zeta) {
            -(nu - 1) * tan(theta) + (1 + nu) * zeta * cos(theta)
        }
        theta.max <- uniroot(fcn, interval = c(-pi/2, pi/2), 
            nu = nu, zeta = zeta)$root
        if (is.finite(lands.cond.t.prop.density.polar(theta.max, 
            nu, zeta))) {
            return(lands.cond.t.prop.density.polar(theta = atan(tau/sqrt(nu)), 
                nu = nu, zeta = zeta)/integrate(lands.cond.t.prop.density.polar, 
                -pi/2, pi/2, nu = nu, zeta = zeta)$integral)
        }
        else {
            quasi.density.fcn.dlands.t <- function(theta, nu, 
                zeta, ln.c) {
                exp(ln.c + (nu - 1) * log(cos(theta)) - (nu/2) * 
                  log(nu) + (1 + nu) * zeta * sin(theta))
            }
            fcn.to.minimize.dlands.t <- function(ln.c, nu, zeta) {
                (integrate(quasi.density.fcn.dlands.t, -pi/2, 
                  pi/2, nu = nu, zeta = zeta, ln.c = ln.c)$value - 
                  1)^2
            }
            start <- -((nu - 1) * log(cos(theta.max)) - (nu/2) * 
                log(nu) + (1 + nu) * zeta * sin(theta.max))
            ln.c <- nlminb(objective = fcn.to.minimize.dlands.t, 
                start = start, nu = nu, zeta = zeta)$par
            return(quasi.density.fcn.dlands.t(theta = atan(tau/sqrt(nu)), 
                nu, zeta, ln.c = ln.c))
        }
    }
}
