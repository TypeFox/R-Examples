num.deriv <- function(f, x, epsilon, ...)
{
    (f(x + epsilon, ...) - f(x - epsilon, ...)) / (2 * epsilon)
}
