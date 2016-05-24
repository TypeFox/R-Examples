Dunnetts.K <-
function (n, df = n - 1, k, m, method = c("Bonferroni", "exact"), 
    type = c("two-sided", "upper", "lower"), conf.level = 0.95) 
{
    method <- match.arg(method)
    type <- match.arg(type)
    if (type == "two-sided") {
        alpha <- (1 - conf.level)/2
    }
    else {
        alpha <- 1 - conf.level
    }
    K <- qt(1 - alpha/k, df) * sqrt(1/m + 1/n)
    if (k > 1 && method == "exact") {
        rho <- 1/((n/m) + 1)
        if (type == "two-sided") {
            fcn.to.min <- function(d, n, df, k, rho, conf.level) {
                (integrate(Dunnetts.K.fcn.to.integrate.2, lower = 0, 
                  upper = Inf, d = d, n = n, df = df, k = k, 
                  rho = rho)$value - conf.level)^2
            }
        }
        else {
            fcn.to.min <- function(d, n, df, k, rho, conf.level) {
                (integrate(Dunnetts.K.fcn.to.integrate.1, lower = 0, 
                  upper = Inf, d = d, n = n, df = df, k = k, 
                  rho = rho)$value - conf.level)^2
            }
        }
        d <- nlminb(start = K/sqrt(1/m + 1/n), objective = fcn.to.min, 
            n = n, df = df, k = k, rho = rho, conf.level = conf.level, 
            lower = .Machine$double.eps)$par
        K <- d * sqrt(1/m + 1/n)
    }
    K
}
