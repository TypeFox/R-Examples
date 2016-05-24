plands.t <-
function (q, nu, zeta, tol = .Machine$double.eps^0.9) 
{
    if (nu < 2) 
        stop("'nu' must be greater than or equal to 2")
    if (q < -.Machine$double.xmax) 
        return(0)
    else if (q > .Machine$double.xmax) 
        return(1)
    else {
        fcn <- function(theta, nu, zeta) {
            -(nu - 1) * tan(theta) + (1 + nu) * zeta * cos(theta)
        }
        theta.max <- uniroot(fcn, interval = c(-pi/2, pi/2), 
            tol = tol, nu = nu, zeta = zeta)$root
        if (is.finite(lands.cond.t.prop.density.polar(theta.max, 
            nu, zeta))) {
            scaling.factor <- integrate(lands.cond.t.prop.density.polar, 
                -pi/2, pi/2, nu = nu, zeta = zeta)$integral
            density.fcn.plands.t <- function(theta, nu, zeta) {
                lands.cond.t.prop.density.polar(theta, nu, zeta)/scaling.factor
            }
        }
        else {
            quasi.density.fcn.plands.t <- function(theta, nu, 
                zeta, ln.c) {
                exp(ln.c + (nu - 1) * log(cos(theta)) - (nu/2) * 
                  log(nu) + (1 + nu) * zeta * sin(theta))
            }
            fcn.to.minimize.plands.t <- function(ln.c, nu, zeta) {
                (integrate(quasi.density.fcn.plands.t, -pi/2, 
                  pi/2, nu = nu, zeta = zeta, ln.c = ln.c)$value - 
                  1)^2
            }
            start <- -((nu - 1) * log(cos(theta.max)) - (nu/2) * 
                log(nu) + (1 + nu) * zeta * sin(theta.max))
            ln.c <- nlminb(objective = fcn.to.minimize.plands.t, 
                start = start, nu = nu, zeta = zeta)$par
            density.fcn.plands.t <- function(theta, nu, zeta) {
                quasi.density.fcn.plands.t(theta, nu, zeta, ln.c = ln.c)
            }
        }
        integrate(density.fcn.plands.t, -pi/2, atan(q/sqrt(nu)), 
            nu = nu, zeta = zeta)$value
    }
}
