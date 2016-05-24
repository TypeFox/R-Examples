# postplnre.lnre.fzm <- function (model, q, m, N, lower.tail=FALSE, ...) 
# {
#   if (! inherits(model, "lnre.fzm")) stop("argument must be object of class 'lnre.fzm'")
#   if (!(is.numeric(N) && all(N >= 0))) stop("argument 'N' must be non-negative integer")
#   if (!(is.numeric(m) && all(m >= 1))) stop("argument 'm' must be positive integer")
#   if (!(is.numeric(q) && all(q >= 0))) stop("argument 'q' must be numeric and non-negative")
# 
#   alpha <- model$param$alpha
#   A <- model$param$A
#   B <- model$param$B
#   C <- model$param2$C
# 
#   rho <- pmin(pmax(q, A), B) # clamp cutoff point rho to model range, so both integrals are valid
# 
#   ## ******* much room for optimisation here *********
# #  termA <- exp(Igamma(m - alpha, N * A, lower=FALSE, log=TRUE) - Cgamma(m + 1, log=TRUE))
#     # termA = Igamma(m - alpha, N * A, lower=FALSE) / m!
# #  termB <- exp(Igamma(m - alpha, N * B, lower=FALSE, log=TRUE) - Cgamma(m + 1, log=TRUE))
#     # termB = Igamma(m - alpha, N * B, lower=FALSE) / m!
# #  termRho <- exp(Igamma(m - alpha, N * rho, lower=FALSE, log=TRUE) - Cgamma(m + 1, log=TRUE))
#     # termRho = Igamma(m - alpha, N * rho, lower=FALSE) / m!
# 
#   ## ****** TODO CHECK: same factor m! in all terms cancels, use Rgamma to avoid overflows
#     termA <- Rgamma(m - alpha, N * A, lower=FALSE)
#     termB <- Rgamma(m - alpha, N * B, lower=FALSE)
#     termRho <- Rgamma(m - alpha, N * rho, lower=FALSE)
# 
#   part <- if (lower.tail) termA - termRho else termRho - termB
#   return(part / (termA - termB)) # factor C * N^alpha / m! should cancel!
# #  denom <- C * N^alpha * (termA - termB)  }
# 
# ## ****** TODO: add similar code for approximate calculation *******
# #  else {
# #    factor <- exp(Igamma(m - alpha, N * A, lower=FALSE, log=TRUE) - Cgamma(m + 1, log=TRUE))
# #    # factor = Igamma(m - alpha, N * A, lower=FALSE) / m!
# #    C * N^alpha * factor
# #  }
# }
