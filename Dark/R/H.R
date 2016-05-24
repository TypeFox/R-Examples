H <- function(x, k = 100, t) {
# x is the measured time
# k is the transition constant, set arbitrarily high
# t is the time at which the transition occurs

	round(1/(1 + exp(-2 * k * (x - t))), 1)
}