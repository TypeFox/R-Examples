# gtlnre <- function (model, m, N=NULL) {
#   if (! inherits(model, "lnre")) stop("first argument must be object of class 'lnre'")
#   if (missing(N)) N <- N(model)
#   if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")
#   if (!(is.numeric(m) && all(m >= 1))) stop("argument 'm' must be positive integer")
# 
#   ## TODO: add code for special case m=0 (with finite population size S, possibly specified by user)
#   (m+1) * EVm(model, m+1, N) / EVm(model, m, N)
# }