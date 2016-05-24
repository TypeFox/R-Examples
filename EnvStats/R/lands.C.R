lands.C <-
function (S, nu, conf.level, ...) 
{
    if (S < .Machine$double.eps) 
        stop("'S' must be positive")
    if (nu < 2) 
        stop("'nu' must be greater than or equal to 2")
    if (conf.level < 0 || conf.level > 1) 
        stop("'conf.level' must be between 0 and 1")
    alpha <- 1 - conf.level
    fcn.to.minimize <- function(m, S, nu, alpha) {
        T.m <- (sqrt(nu + 1) * ((-S^2)/2 - m))/S
        zeta.m <- (-S * sqrt(nu + (T.m^2)))/(2 * sqrt(nu + 1))
        (T.m - qlands.t(p = alpha, nu = nu, zeta = zeta.m))^2
    }
    start.T.m <- (-sqrt(nu + 1) * S)/2
    start.zeta.m <- (-S * sqrt(nu + (start.T.m^2)))/(2 * sqrt(nu + 
        1))
    start.T.m <- qlands.t(1 - conf.level, nu, start.zeta.m)
    start.m <- ((-S^2)/2) - (S * start.T.m)/sqrt(nu + 1)
    m <- nlminb(start = start.m, objective = fcn.to.minimize, 
        ..., S = S, nu = nu, alpha = alpha)$par
    (m * sqrt(nu))/S
}
