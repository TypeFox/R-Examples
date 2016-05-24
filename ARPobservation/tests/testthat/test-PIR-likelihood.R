context("Evaluating the functioning of the likelihood estimator")

library(plyr)
library(reshape2)

test_that("As phi gets very close to zero, the psi values become constant and the log likelihood approaches a known value",{
  
  phi <- 10^seq(-3,-6,-1)
  zeta <- seq(.10,1, length.out = 9)
  set.seed(1028194800)
  U <- sample(0:1, 30, replace = TRUE)
  c <- 1
  d <- 5
  
  parms <- expand.grid(phi = phi, zeta = zeta)
  
  results <- mdply(.data = parms, .fun = PIR_loglik, U = U, c = c, d = d, .inform = TRUE)
  results_wide <- dcast(results, zeta ~ -phi, value.var = "V1")
  
  
  theoretical <- ((log(1-exp(-1 * zeta * c)) + zeta * c)* sum(U)) - (length(U) * zeta * c)
  
  all_values <- cbind(results_wide, theoretical)
  
  psi <- mdply(.data = phi, .fun = PIRpsi, zeta = .25, U = U, c = c, d=d)
  
  expect_true(all(apply(psi[,-1], 1, function(x) max(x) - min(x)) < sqrt(.Machine$double.eps)))
  expect_that(all(abs(all_values[,2] - all_values[,3]) > abs(all_values[,3] - all_values[,4]) & 
                    abs(all_values[,3] - all_values[,4]) > abs(all_values[,4] - all_values[,5]) & 
                    abs(all_values[,4] - all_values[,5]) > abs(all_values[,5] - all_values[,6]) & 
                    abs(all_values[,5] - all_values[,6]) < 10^-3), is_true())  
})

test_that("If d is much larger than 1/zeta, then psis = phi and the log-likelihood approaches a known value",{
  phi <- .25
  zeta <- seq(.10,1, length.out = 9)
  set.seed(1028194800)
  U <- sample(0:1, 30, replace = TRUE)
  c <- 1
  d <- c(2 / zeta, 3 / zeta, 4 / zeta, 5 / zeta, 10 / zeta)
  
  parms <- cbind(zeta, d)
  
  results <- mdply(.data = parms, .fun = PIR_loglik, phi = phi, U = U, c = c, .inform = TRUE)
  results$d <- with(results, d*zeta)
  results_wide <- dcast(results, zeta ~ d, value.var = "V1")
  
  theoretical <- (log(1 / (1 - phi) - exp(-1 * zeta * c/(1-phi))) + zeta * c / (1-phi)) * sum(U) + length(U) * (log(1 - phi) - zeta * c / (1-phi))
  
  all_values <- cbind(results_wide, theoretical)
  psi <- mdply(.data = parms[37:45,], .fun = PIRpsi, phi = .25, U = U, c = c)
  
  expect_true(all(apply(psi[,-2:-1], 1, function(x) max(x) - min(x)) < sqrt(.Machine$double.eps)))
  
  expect_that(all(abs(all_values[,2] - all_values[,3]) > abs(all_values[,3] - all_values[,4]) & 
                    abs(all_values[,3] - all_values[,4]) > abs(all_values[,4] - all_values[,5]) & 
                    abs(all_values[,4] - all_values[,5]) > abs(all_values[,5] - all_values[,6]) & 
                    abs(all_values[,5] - all_values[,6]) < 10^-3), is_true())
})

test_that("If all Us = 0, then the loglikelihood is some known value", {
  phi <- .75
  zeta <- seq(.10,1, length.out = 9)
  U <- rep(0,30)
  c <- 1
  d <- 0
  
  loglik <-ldply(.data = zeta, .fun = PIR_loglik, phi = phi, U = U, c = c, d = d)
  
  theoretical <- log(1-phi) - zeta * c /(1-phi) + (length(U) - 1) * (log(1 - p_0(d, phi, zeta)) - zeta * c/(1-phi))
  
  expect_that(loglik$V1, equals(theoretical))
})