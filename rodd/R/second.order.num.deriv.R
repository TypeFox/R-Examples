second.order.num.deriv <- function(f, x, epsilon, ...)
{
	(f(x + epsilon, ...) - 2 * f(x, ...) + f(x - epsilon, ...)) / (epsilon * epsilon)
}

