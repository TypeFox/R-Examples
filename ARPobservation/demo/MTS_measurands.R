set.seed(10)
# set parameters
zeta <- 1 / 60
iterations <- 100
phi <- rep(seq(0.01, 0.99, 0.01), each = iterations)
mu <- phi / zeta 
lambda <- (1 - phi) / zeta
# generate behavior streams, apply continuous recording and MTS
BS <- r_behavior_stream(n = length(phi), mu = mu, lambda = lambda, 
                        F_event = F_exp(), F_interim = F_exp(), stream_length = 600)
obs <- reported_observations(BS, c("C","M"), interval_length = 15)
# plot results
library(ggplot2)
qplot(x = phi, y = M, data = obs, geom = "point", alpha = I(0.1),
      xlab = "Prevalence", ylab = "Momentary time sampling") + 
  geom_smooth(method = "loess", se = FALSE, color = "red") + 
  theme_bw()
qplot(x = C, y = M, data = obs, geom = "point", alpha = I(0.1),
      xlab = "Continuous recording", ylab = "Momentary time sampling") +
  geom_smooth(method = "loess", se = FALSE, color = "red") + 
  theme_bw()
