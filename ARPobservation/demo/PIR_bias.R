set.seed(10)
# set parameters
stream_length <- 60
iterations <- 1000
phi <- 0.1
zeta <- seq(1/100, 1/2, length.out = 100)
interim_dist <- c("Weibull(3)", "Gamma(3)", "Uniform",
                  "Exponential", "Gamma(1/2)", "Weibull(1/2)")
dist_lookup <- list("Weibull(3)" = F_weib(3), "Gamma(3)" = F_gam(3),
                    "Uniform" = F_unif(), "Exponential" = F_exp(),
                    "Gamma(1/2)" = F_gam(1/2), "Weibull(1/2)" = F_weib(1/2))
# create parameter combinations
parms <- expand.grid(interim_dist = interim_dist, zeta = zeta)
# function for simulating mean of PIR data
PIR_mean <- function(phi, zeta, interim_dist, iterations, stream_length) {
  BS <- r_behavior_stream(n = iterations, 
                          mu = phi / zeta, lambda = (1 - phi) / zeta, 
                          F_event = F_const(), 
                          F_interim = dist_lookup[[interim_dist]], 
                          stream_length = stream_length)
  PIR <- interval_recording(BS, interval_length = 1)
  c(PIR_mean = mean(PIR))
}
# apply function to each combination of parameter values
library(plyr)
PIR_sim <- mdply(parms, PIR_mean, phi=phi, 
                 iterations = iterations, stream_length = stream_length)
# plot results
library(ggplot2)
qplot(x = zeta, y = PIR_mean, color = interim_dist, linetype = interim_dist,
      data = PIR_sim, geom = "smooth", method = "loess", se = FALSE,
      xlab = "Incidence (per interval)", ylab = expression(E(Y^P))) + 
  labs(linetype="Interim time distribution", color="Interim time distribution") +
  theme_bw()
