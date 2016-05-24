"fun.percentile" <-
function(x, p)
{
if(p < 0 || p > 1) {
stop("p must be between 0 and 1")
}
n <- length(x)
x <- sort(x)
if(p < 1/(n + 1) || p > n/(n + 1)) {
stop("p does not exist")
}
E <- (n + 1) * p
r <- trunc(E)
ab <- E - r
x[r] + ab * (x[r + 1] - x[r])
}

