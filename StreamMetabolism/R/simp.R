`simp` <- function (y, a = NULL, b = NULL, x = NULL, n = 200)
{
   if (is.null(a) | is.null(b)) {
       if (is.null(x))
           stop("No x values provided to integrate over.\n")
   }
   else {
       x <- c(a, b)
   }
   fff <- 1
   if (length(x) == 2) {
       if (x[1] == x[2])
           return(0)
       if (x[2] < x[1]) {
           fff <- -1
           x <- rev(x)
       }
       x <- seq(x[1], x[2], length = n)
       if (is.function(y))
           y <- y(x)
       else {
           cat("y must be a function when x is\n")
           cat("of length equal to 2.\n")
           stop("Bailing out.\n")
       }
       equisp <- TRUE
   }
   else {
       if (is.function(y))
           y <- y(x)
       else if (length(y) != length(x))
           stop("Mismatch in lengths of x and y.\n")
       s <- order(x)
       x <- x[s]
       ddd <- diff(x)
       if (any(ddd == 0))
           stop("Gridpoints must be distinct.\n")
       equisp <- isTRUE(all.equal(diff(ddd), rep(0, length(ddd) - 1)))
       y <- y[s]
   }
   n <- length(x) - 1
   if (equisp) {
       old.op <- options(warn = -1)
       on.exit(options(old.op))
       M <- matrix(y, nrow = n + 2, ncol = 4)[1:(n - 2), ]
       h <- x[2] - x[1]
       fc <- h * c(-1, 13, 13, -1)/24
       aa <- apply(t(M) * fc, 2, sum)
       a1 <- h * sum(y[1:3] * c(5, 8, -1))/12
       an <- h * sum(y[(n - 1):(n + 1)] * c(-1, 8, 5))/12
       return(fff * sum(c(a1, aa, an)))
   }
   m <- n%/%2
   i <- 1:(m + 1)
   a <- x[2 * i] - x[2 * i - 1]
   i <- 1:m
   b <- x[2 * i + 1] - x[2 * i]
   o <- (a[i] * b + 2 * a[i] * a[i] - b * b)/(6 * a[i])
   p <- (a[i] + b)^3/(6 * a[i] * b)
   q <- (a[i] * b + 2 * b * b - a[i] * a[i])/(6 * b)
   k <- numeric(n + 1)
   k[1] <- o[1]
   i <- 1:(m - 1)
   k[2 * i] <- p[i]
   k[2 * i + 1] <- q[i] + o[-1]
   if (n > 2 * m) {
       aa <- a[m + 1]
       bb <- b[m]
       den <- 6 * bb * (bb + aa)
       k[2 * m] <- p[m] - (aa^3)/den
       k[2 * m + 1] <- q[m] + (aa^3 + 4 * bb * aa^2 + 3 * aa *
           bb^2)/den
       k[2 * m + 2] <- (2 * bb * aa^2 + 3 * aa * bb^2)/den
   }
   else {
       k[2 * m] <- p[m]
       k[2 * m + 1] <- q[m]
   }
   fff * sum(k * y)
}