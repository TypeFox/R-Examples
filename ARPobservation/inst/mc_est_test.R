#this unit test is "optional" because the small simulation contained within it would significantly lengthen the build time every time the package was checked

library(testthat)
library(ARPobservation)
library(plyr)

Context("An optional unit test for examining the log-likelihood estimator.")


test_that("Monte Carlo estimates match the expected value of the likelihood function", {

phi <- seq(.10, .50, length.out = 3)
zeta <- seq(.1, .50, length.out = 3)
c <- 1
d <- c(0, 1/4, 1/3)
iterations <- 100000
nIntervals <- 5
set.seed(847291)

parms <- expand.grid(phi = phi, zeta = zeta, c = c, d = d, 
                     nIntervals = nIntervals, iterations = iterations)

generateEsts <- function(phi, zeta, c, d, nIntervals, iterations) {  
  mu <- phi/zeta
  lambda <- (1-phi)/zeta
  intervalLength <- c + d
  
  BS <- r_behavior_stream(n = iterations, mu = mu, lambda = lambda, 
                          F_event = F_exp(), F_interim = F_exp(), 
                          stream_length = intervalLength * nIntervals)
  samples <- interval_recording(BS = BS, interval_length = intervalLength, rest_length = d, summarize = FALSE)
  
  samples <- as.data.frame(t(samples))
  samples$One <- 1
  mc_est <- ddply(samples, .variables = names(samples)[-(nIntervals + 1)], summarize, freq = sum(One)) 
  mc_est$p <- exp(apply(mc_est[,1:nIntervals], 1, PIR_loglik, phi = phi, zeta = zeta, c = c, d = d))
  return (mc_est)
}

mc_est <- mdply(.data = parms, .fun = generateEsts, .inform = TRUE)

mc_est<- mc_est[,c(1,2,4,12,13)]

chisq_wrapper <- function(subset){
 chisq.test(x = subset$freq, p = subset$p)$p.value
}
test_values <- ddply(mc_est, .(phi, zeta, d), .fun = chisq_wrapper)

expect_that(min(test_values$V1) > .01, is_true())
})
