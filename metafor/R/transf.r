transf.abt <-
function (xi, ...) 
{
    zi <- -log(1 - xi)
    return(c(zi))
}
transf.ahw <-
function (xi, ...) 
{
    zi <- 1 - (1 - xi)^(1/3)
    return(c(zi))
}
transf.arcsin <-
function (xi, ...) 
{
    zi <- asin(sqrt(xi))
    return(c(zi))
}
transf.exp.int <-
function (xi, targs = NULL, ...) 
{
    if (is.null(targs$tau2)) 
        targs$tau2 <- 0
    if (is.null(targs$lower)) 
        targs$lower <- xi - 5 * sqrt(targs$tau2)
    if (is.null(targs$upper)) 
        targs$upper <- xi + 5 * sqrt(targs$tau2)
    toint <- function(zval, xi, tau2) {
        exp(zval) * dnorm(zval, mean = xi, sd = sqrt(tau2))
    }
    cfunc <- function(xi, tau2, lower, upper) {
        integrate(toint, lower = lower, upper = upper, xi = xi, 
            tau2 = tau2)$value
    }
    zi <- mapply(xi, FUN = cfunc, tau2 = targs$tau2, lower = targs$lower, 
        upper = targs$upper)
    return(c(zi))
}
transf.iabt <-
function (xi, ...) 
{
    zi <- 1 - exp(-xi)
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi < 0] <- 0
    return(c(zi))
}
transf.iahw <-
function (xi, ...) 
{
    zi <- 1 - (1 - xi)^3
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi > 1] <- 1
    zi[xi < 0] <- 0
    return(c(zi))
}
transf.iarcsin <-
function (xi, ...) 
{
    zi <- sin(xi)^2
    zi[xi < 0] <- 0
    zi[xi > asin(1)] <- 1
    return(c(zi))
}
transf.iirft <-
function (xi, ti, ...) 
{
    zi <- (1/ti - 8 * xi^2 + 16 * ti * xi^4)/(16 * xi^2 * ti)
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi < transf.irft(0, ti)] <- 0
    zi[zi <= .Machine$double.eps] <- 0
    return(c(zi))
}
transf.ilogit <-
function (xi, ...) 
{
    zi <- exp(xi)/(1 + exp(xi))
    zi[xi == -Inf] <- 0
    zi[xi == Inf] <- 1
    zi[is.nan(zi) & (xi > 0.5)] <- 1
    return(c(zi))
}
transf.ilogit.int <-
function (xi, targs = NULL, ...) 
{
    if (is.null(targs$tau2)) 
        targs$tau2 <- 0
    if (is.null(targs$lower)) 
        targs$lower <- xi - 5 * sqrt(targs$tau2)
    if (is.null(targs$upper)) 
        targs$upper <- xi + 5 * sqrt(targs$tau2)
    toint <- function(zval, xi, tau2) {
        zi <- exp(zval)/(1 + exp(zval))
        zi[xi == -Inf] <- 0
        zi[xi == Inf] <- 1
        zi[is.nan(zi) & (xi > 0.5)] <- 1
        zi * dnorm(zval, mean = xi, sd = sqrt(tau2))
    }
    cfunc <- function(xi, tau2, lower, upper) {
        integrate(toint, lower = lower, upper = upper, xi = xi, 
            tau2 = tau2)$value
    }
    zi <- mapply(xi, FUN = cfunc, tau2 = targs$tau2, lower = targs$lower, 
        upper = targs$upper)
    return(c(zi))
}
transf.ipft <-
function (xi, ni, ...) 
{
    zi <- 1/2 * (1 - sign(cos(2 * xi)) * sqrt(1 - (sin(2 * xi) + 
        (sin(2 * xi) - 1/sin(2 * xi))/ni)^2))
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi > transf.pft(1, ni)] <- 1
    zi[xi < transf.pft(0, ni)] <- 0
    return(c(zi))
}
transf.ipft.hm <-
function (xi, targs, ...) 
{
    o.warn <- getOption("warn")
    on.exit(options(warn = o.warn))
    options(warn = -1)
    ni <- 1/(mean(1/targs$ni, na.rm = TRUE))
    zi <- 1/2 * (1 - sign(cos(2 * xi)) * sqrt(1 - (sin(2 * xi) + 
        (sin(2 * xi) - 1/sin(2 * xi))/ni)^2))
    zi <- ifelse(is.nan(zi), NA, zi)
    zi[xi > transf.pft(1, ni)] <- 1
    zi[xi < transf.pft(0, ni)] <- 0
    return(c(zi))
}
transf.irft <-
function (xi, ti, ...) 
{
    zi <- 1/2 * (sqrt(xi) + sqrt(xi + 1/ti))
    return(c(zi))
}
transf.isqrt <-
function (xi, ...) 
{
    zi <- xi * xi
    zi[xi < 0] <- 0
    return(c(zi))
}
transf.logit <-
function (xi, ...) 
{
    zi <- log(xi/(1 - xi))
    return(c(zi))
}
transf.pft <-
function (xi, ni, ...) 
{
    xi <- xi * ni
    zi <- 1/2 * (asin(sqrt(xi/(ni + 1))) + asin(sqrt((xi + 1)/(ni + 
        1))))
    return(c(zi))
}
transf.rtoz <-
function (xi, ...) 
{
    zi <- 1/2 * log((1 + xi)/(1 - xi))
    return(c(zi))
}
transf.ztor <-
function (xi, ...) 
{
    zi <- (exp(2 * xi) - 1)/(exp(2 * xi) + 1)
    zi[xi == -Inf] <- -1
    zi[xi == Inf] <- 1
    zi[is.nan(zi) & (xi > 0)] <- 1
    return(c(zi))
}
transf.ztor.int <-
function (xi, targs = NULL, ...) 
{
    if (is.null(targs$tau2)) 
        targs$tau2 <- 0
    if (is.null(targs$lower)) 
        targs$lower <- xi - 5 * sqrt(targs$tau2)
    if (is.null(targs$upper)) 
        targs$upper <- xi + 5 * sqrt(targs$tau2)
    toint <- function(zval, xi, tau2) {
        zi <- (exp(2 * zval) - 1)/(exp(2 * zval) + 1)
        zi[xi == -Inf] <- -1
        zi[xi == Inf] <- 1
        zi[is.nan(zi) & (xi > 0)] <- 1
        zi * dnorm(zval, mean = xi, sd = sqrt(tau2))
    }
    cfunc <- function(xi, tau2, lower, upper) {
        integrate(toint, lower = lower, upper = upper, xi = xi, 
            tau2 = tau2)$value
    }
    zi <- mapply(xi, FUN = cfunc, tau2 = targs$tau2, lower = targs$lower, 
        upper = targs$upper)
    return(c(zi))
}
