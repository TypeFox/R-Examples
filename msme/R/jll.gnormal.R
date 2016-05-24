
jll_gnormal <- function(params, y, X, group, ...) {
   p <- ncol(X)
   N_i <- tapply(y, group, length)
   Su <- exp(params[p+1])
   Se <- exp(params[p+2])
   z <- y - X %*% params[1:p]
   gamma_i <- Su^2 / (N_i * Su^2 + Se^2)
   c1 <- (tapply(z^2, group, sum) -
          gamma_i * tapply(z, group, sum)^2) / Se^2
   c2 <- log(N_i * Su^2 / Se^2 + 1)
   c3 <- N_i * log(2 * pi * Se^2)
   return(sum(-0.5 * (c1 + c2 + c3)))
}
