
 # test batching (blen) and spacing (nspac) together

 epsilon <- 1e-15

 library(mcmc)

 RNGkind("Marsaglia-Multicarry")
 set.seed(42)

 n <- 100
 rho <- 0.5
 beta0 <- 0.25
 beta1 <- 1
 beta2 <- 0.5

 x1 <- rnorm(n)
 x2 <- rho * x1 + sqrt(1 - rho^2) * rnorm(n)
 eta <- beta0 + beta1 * x1 + beta2 * x2
 p <- 1 / (1 + exp(- eta))
 y <- as.numeric(runif(n) < p)

 out <- glm(y ~ x1 + x2, family = binomial())

 logl <- function(beta) {
     if (length(beta) != 3) stop("length(beta) != 3")
     beta0 <- beta[1]
     beta1 <- beta[2]
     beta2 <- beta[3]
     eta <- beta0 + beta1 * x1 + beta2 * x2
     p <- exp(eta) / (1 + exp(eta))
     return(sum(log(p[y == 1])) + sum(log(1 - p[y == 0])))
 }

 out.metro <- metrop(logl, coefficients(out), 1e3, scale = 0.01)
 out.metro$accept

 out.metro <- metrop(out.metro, scale = 0.1)
 out.metro$accept

 out.metro <- metrop(out.metro, scale = 0.5)
 out.metro$accept

 apply(out.metro$batch, 2, mean)

 out.metro <- metrop(logl, as.numeric(coefficients(out)), 1e2,
     scale = 0.5, debug = TRUE, blen = 5, nspac = 3)

 niter <- out.metro$nbatch * out.metro$blen * out.metro$nspac
 niter == nrow(out.metro$current)
 niter == nrow(out.metro$proposal)
 all(out.metro$current[1, ] == out.metro$initial)
 all(out.metro$current[niter, ] == out.metro$final) |
     all(out.metro$proposal[niter, ] == out.metro$final)

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

 my.curr.log.green <- apply(out.metro$current, 1, logl)
 my.prop.log.green <- apply(out.metro$proposal, 1, logl)
 all(is.na(out.metro$u) == (my.prop.log.green > my.curr.log.green))
 foo <- my.prop.log.green - my.curr.log.green
 max(abs(foo - out.metro$log.green)) < epsilon

 my.accept <- is.na(my.u) | my.u < exp(foo)
 sum(my.accept) == round(n * out.metro$accept)
 if (my.accept[niter]) {
     all(out.metro$proposal[niter, ] == out.metro$final)
 } else {
     all(out.metro$current[niter, ] == out.metro$final)
 }

 my.current <- out.metro$current
 my.current[my.accept, ] <- my.proposal[my.accept, ]
 my.current <- rbind(out.metro$initial, my.current[- niter, ])
 max(abs(out.metro$current - my.current)) < epsilon

 my.path <- matrix(NA, n, d)
 my.path[my.accept, ] <- out.metro$proposal[my.accept, ]
 my.path[! my.accept, ] <- out.metro$current[! my.accept, ]
 nspac <- out.metro$nspac

 my.path <- my.path[seq(nspac, niter, by = nspac), ]

 foom <- array(as.vector(t(my.path)), c(d, out.metro$blen, out.metro$nbatch))
 boom <- t(apply(foom, c(1, 3), mean))

 all(dim(boom) == dim(out.metro$batch))
 max(abs(boom - out.metro$batch)) < epsilon

