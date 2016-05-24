f_rastrigin <-
function(x) {
d <- length(x)
10*d+sum(x*x-10*cos(2*pi*x))
}
