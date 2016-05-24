qlands.t <-
function (p, nu, zeta, tol = .Machine$double.eps^0.9) 
{
    if (p < 0 || p > 1) 
        stop("'p' must be between 0 and 1")
    if (nu < 2) 
        stop("'nu' must be greater than or equal to 2")
    if (p == 0) 
        return(-Inf)
    else if (p == 1) 
        return(Inf)
    else {
        fcn <- function(theta, nu, zeta) {
            -(nu - 1) * tan(theta) + (1 + nu) * zeta * cos(theta)
        }
        theta.max <- uniroot(fcn, interval = c(-pi/2, pi/2), 
            tol = tol, nu = nu, zeta = zeta)$root
        if (is.finite(lands.cond.t.prop.density.polar(theta.max, 
            nu, zeta))) {
            scaling.factor <- integrate(lands.cond.t.prop.density.polar, 
                -pi/2, pi/2, nu = nu, zeta = zeta)$value
            density.fcn.qlands.t <- function(theta, nu, zeta) {
                lands.cond.t.prop.density.polar(theta, nu, zeta)/scaling.factor
            }
        }
        else {
            quasi.density.fcn.qlands.t <- function(theta, nu, 
                zeta, ln.c) {
                exp(ln.c + (nu - 1) * log(cos(theta)) - (nu/2) * 
                  log(nu) + (1 + nu) * zeta * sin(theta))
            }
            fcn.to.minimize.qlands.t <- function(ln.c, nu, zeta) {
                (integrate(quasi.density.fcn.qlands.t, -pi/2, 
                  pi/2, nu = nu, zeta = zeta, ln.c = ln.c)$value - 
                  1)^2
            }
            start <- -((nu - 1) * log(cos(theta.max)) - (nu/2) * 
                log(nu) + (1 + nu) * zeta * sin(theta.max))
            ln.c <- nlminb(objective = fcn.to.minimize.qlands.t, 
                start = start, nu = nu, zeta = zeta)$par
            density.fcn.qlands.t <- function(theta, nu, zeta) {
                quasi.density.fcn.qlands.t(theta, nu, zeta, ln.c = ln.c)
            }
        }
        uniroot.fcn.qlands.t <- function(theta, nu, zeta, p) {
            integrate(density.fcn.qlands.t, -pi/2, theta, nu = nu, 
                zeta = zeta)$value - p
        }
        if (p > 0.5) {
            new.zeta <- -zeta
            new.p <- 1 - p
        }
        else {
            new.zeta <- zeta
            new.p <- p
        }
        quant <- uniroot(uniroot.fcn.qlands.t, interval = c(-pi/2, 
            pi/2), tol = tol, nu = nu, zeta = new.zeta, p = new.p)$root
        quant <- sqrt(nu) * tan(quant)
        if (p > 0.5) 
            quant <- -quant
        quant
    }
}
