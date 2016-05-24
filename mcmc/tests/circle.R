
 epsilon <- 1e-15

 library(mcmc)

 RNGkind("Marsaglia-Multicarry")
 set.seed(42)

 d <- 5

 logh <- function(x) {
     if (! is.numeric(x)) stop("x not numeric")
     if (length(x) != d) stop("length(x) != d")
     fred <- 1 - sum(x^2)
     if (fred > 0) return(log(fred)) else return(-Inf)
 }

 out.metro <- metrop(logh, rep(0, d), 1e3, scale = 0.01)
 out.metro$accept

 out.metro <- metrop(out.metro, scale = 0.1)
 out.metro$accept

 out.metro <- metrop(out.metro, scale = 0.5)
 out.metro$accept

 out.metro <- metrop(out.metro, scale = 0.4)
 out.metro$accept

 out.metro <- metrop(out.metro, nbatch = 1e2, debug = TRUE)

 all(out.metro$batch[- out.metro$nbatch, ] == out.metro$current[- 1, ])
 all(out.metro$current[1, ] == out.metro$initial)
 all(out.metro$batch[out.metro$nbatch, ] == out.metro$final)

 .Random.seed <- out.metro$initial.seed
 d <- ncol(out.metro$proposal)
 n <- nrow(out.metro$proposal)
 my.proposal <- matrix(NA, n, d)
 my.u <- double(n)
 ska <- out.metro$scale
 for (i in 1:n) {
     my.proposal[i, ] <- out.metro$current[i, ] + ska * rnorm(d)
     if (is.na(out.metro$u[i])) {
         my.u[i] <- NA
     } else {
         my.u[i] <- runif(1)
     }
 }
 max(abs(out.metro$proposal - my.proposal)) < epsilon
 all(is.na(out.metro$u) == is.na(my.u))
 all(out.metro$u[!is.na(out.metro$u)] == my.u[!is.na(my.u)])

 my.curr.log.green <- apply(out.metro$current, 1, logh)
 my.prop.log.green <- apply(out.metro$proposal, 1, logh)
 all(is.na(out.metro$u) == ((my.prop.log.green == -Inf) |
     (my.prop.log.green > my.curr.log.green)))
 foo <- my.prop.log.green - my.curr.log.green
 blurfle <- foo - out.metro$log.green
 blurfle[foo == -Inf & out.metro$log.green == -Inf] <- 0
 max(blurfle) < epsilon

 my.accept <- (my.prop.log.green > -Inf) & (is.na(my.u) | my.u < exp(foo))
 sum(my.accept) == round(n * out.metro$accept)

 my.path <- matrix(NA, n, d)
 my.path[my.accept, ] <- out.metro$proposal[my.accept, ]
 my.path[! my.accept, ] <- out.metro$current[! my.accept, ]

 all(my.path == out.metro$batch)

