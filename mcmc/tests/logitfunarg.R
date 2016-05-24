
 # test outfun (function)

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
     scale = 0.5, debug = TRUE, outfun = function(x) c(x, x^2))

 out.metro <- metrop(out.metro)
 out.metro$outfun
 dim(out.metro$batch)

 logl <- function(beta, x1, x2, y) {
     if (length(beta) != 3) stop("length(beta) != 3")
     beta0 <- beta[1]
     beta1 <- beta[2]
     beta2 <- beta[3]
     eta <- beta0 + beta1 * x1 + beta2 * x2
     p <- exp(eta) / (1 + exp(eta))
     return(sum(log(p[y == 1])) + sum(log(1 - p[y == 0])))
 }

 out.metro <- metrop(logl, as.numeric(coefficients(out)), 1e2,
     scale = 0.5, debug = TRUE, x1 = x1, x2 = x2, y = y)
 out.metro$lud
 out.metro <- metrop(out.metro, x1 = x1, x2 = x2, y = y)
 out.metro$lud

