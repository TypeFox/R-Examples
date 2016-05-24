## 
## unit test provided by Olivia Lau
## 
cat("Checking cpolr against monte carlo simulations for vector responses\n")
data(freedom)
C <- anchors(self ~ vign1 + vign3 + vign6, data = freedom, method="C")
freedom <- insert(freedom, C)

cat("Generating some sensible parameter starting values\n")
minitrue <- cpolr(cbind(Cs, Ce) ~ as.factor(country) + sex + educ, data = freedom)
beta <- minitrue$coef
tau <- minitrue$zeta
check3 <- c(beta, tau)
x <- model.matrix(minitrue)[, -1]
mu <- x %*% beta

cat("Function to do one monte carlo simulation of model parameters
     tau     = rescaled cutpoints (on latent scale)
     mu      = means for latent variable
     x       = model.matrix
     prob    = probability of censoring
     lambda  = parameter for Poisson distribution
               (for range of vector responses)\n")
montecarlo <- function(tau, mu, x, prob, lambda) {
  ## generating scalar response
  tmp1 <- rnorm(length(mu), mean = mu, sd = 1)
  out <- matrix(FALSE, nrow = length(mu), ncol = length(tau))
  for (i in 1:length(tau)) out[,i] <- tmp1 >= tau[i]
  y <- apply(out, 1, sum)
  Cs <- Ce <- y + 1
  names(Cs) <- names(Ce) <- rownames(x)
  ## idx = indicator for which observaitons are censored
  idx <- NULL
  cidx <- sort(unique(Cs))
  for (i in cidx[-length(cidx)]) {
    tmp2 <- rownames(x)[which(Cs == i)]
    idx <- c(idx, sample(tmp2, size = floor(prob * length(tmp2))))
  }
  ## censoring some scalar responses to ranges
  Cs[idx] <- pmax(Cs[idx] - rpois(length(idx), lambda=lambda), 1)
  Ce[idx] <- pmin(Ce[idx] + rpois(length(idx), lambda=lambda), length(tau) + 1)
  dta <- data.frame(Cs=Cs, Ce=Ce, x = x)
  ## fitting model, returning simulations
  out <- cpolr(cbind(Cs, Ce) ~ x, data = dta)
  c(out$coef, out$zeta)
}

sims <- as.integer(readline
                   ("Select number of monte carlo simulations <then press enter>: "))
if (is.na(sims)) sims <- 50

cat("Running monte carlo simultions\n")
test3 <- replicate(sims, montecarlo(tau = tau, mu = mu, x = x,
                                    prob = 0.1, lambda=0.6))

cat("Matrix with the truth in the first column, the estimated
paramters in the 2nd, and the SD of the simulations in the 3rd.\n")
cbind(check3, apply(test3, 1, mean), apply(test3, 1, sd))




