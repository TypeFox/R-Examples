gaussian.hypergeometric <-
function (a, b, c, x) 
{
    if (!is.vector(a, mode = "numeric") || length(a) != 1 || 
        !is.vector(b, mode = "numeric") || length(b) != 1 || 
        !is.vector(c, mode = "numeric") || length(c) != 1 || 
        !is.vector(x, mode = "numeric") || length(x) != 1) 
        stop("'a', 'b', 'c', and 'x' must be numeric scalars")
    if (a < 0 || c < a) 
        stop("'a' and 'c' must satisfy c > a > 0")
    if (abs(x) >= 1) 
        stop("'x' must be a scalar between -1 and 1")
    if (x == 0) 
        ret.val <- 1
    else if (0 < x && x < 1 && c == a + 1 && b < 1) 
        ret.val <- (a/x^a) * pbeta(x, a, 1 - b)
    else {
        integrand <- function(u, a, b, c, x) {
            u^(a - 1) * (1 - u)^(c - a - 1) * (1 - x * u)^(-b)
        }
        ret.val <- exp(lgamma(c) - lgamma(a) - lgamma(c - a)) * 
            integrate(integrand, lower = 0, upper = 1, a = a, 
                b = b, c = c, x = x)$value
    }
    ret.val
}
