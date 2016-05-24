.normalDensity2d <- function(mu1 = 150, mu2 = 165, sigma1 = 1, sigma2 = 1, rho = 0, x = seq(from = -1, to = +1, length = 10), y = seq(from = -1, 
    to = +1, length = 10)) {
    z = matrix(nrow = length(x), ncol = length(y))
    zx = 0
    zy = 0
    for (x1 in x) {
        zx = zx + 1
        for (x2 in y) {
            zy = zy + 1
            z[zx, zy] = 1/(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) * exp(-1/(2 * (1 - rho^2)) * (((x1 - mu1)/sigma1)^2 - 2 * rho * (x1 - mu1)/sigma1 * (x2 - 
                mu2)/sigma2 + ((x2 - mu2)/sigma2)^2))
        }
        zy = 0
    }
    return(z)
}
simProc = function(x1, x2, x3, noise = TRUE) {
    .norm2d <- function(x1, x2, mu1 = 160, mu2 = 165, rho = 0.7, sigma1 = 45, sigma2 = 22.5) {
        z = 1/(2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) * exp(-1/(2 * (1 - rho^2)) * (((x1 - mu1)/sigma1)^2 - 2 * rho * (x1 - mu1)/sigma1 * (x2 - mu2)/sigma2 + 
            ((x2 - mu2)/sigma2)^2))
        return(z)
    }
    max_z = 0.0002200907
    min_z = 8.358082e-10
    yield = .norm2d(x1 = x1, x2 = x2)
    yield = yield - min_z
    yield = (yield/max_z) * 0.9
    if (noise) 
        yield = yield + rnorm(length(yield), mean = 0, sd = 0.007)
    return(yield)
}
print.invisible <- function(x, ...) {
    cat("")
}
class(simProc) = "invisible" 
