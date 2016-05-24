# postdlnre.lnre <- function (model, x, m, N, ...) {
#   if (! inherits(model, "lnre")) stop("first argument must be object of class 'lnre'")
#   if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")
#   if (!(is.numeric(m) && all(m >= 1))) stop("argument 'm' must be positive integer")
# 
#   factor <- exp( m * log(N * x) - N * x - Cgamma(m + 1, log=TRUE) ) # = (Nx)^m * exp(-Nx) / m!
#   factor * tdlnre(model, x) / EVm(model, m, N) ## ******* was dlnre(), but shouldn't this be tdlnre() instead?? *******
# }
# 
# postldlnre.lnre <- function (model, x, m, N, base=10, log.x=FALSE, ...)
# {
#   if (! inherits(model, "lnre")) stop("first argument must be object of class 'lnre'")
# 
#   if (log.x) x <- base ^ x
#   log(base) * x * postdlnre(model, x, m, N, ...)
# }

