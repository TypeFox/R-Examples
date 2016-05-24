set.seed(10)
# set parameters
zeta <- 1 / 60
iterations <- 5000
phi <- seq(0.01, 0.99, 0.01)
gen_dist <- c("Exponential-Exponential","Gamma(3)-Gamma(3)","Constant-Gamma(3)")
dist_lookup <- list("Exponential-Exponential" = list(F_exp(),F_exp()),
                    "Gamma(3)-Gamma(3)" = list(F_gam(3), F_gam(3)),
                    "Constant-Gamma(3)" = list(F_const(), F_gam(3)))
# create parameter combinations
parms <- expand.grid(gen_dist = gen_dist, phi = phi)
# function for simulating mean and variance of continuous recording data
CDR_moments <- function(phi, zeta, gen_dist, iterations, stream_length) {
  BS <- r_behavior_stream(n = iterations, 
                          mu = phi / zeta, lambda = (1 - phi) / zeta, 
                          F_event = dist_lookup[[gen_dist]][[1]], 
                          F_interim = dist_lookup[[gen_dist]][[2]], 
                          stream_length = stream_length)
  CDR <- continuous_duration_recording(BS)
  c(mean=mean(CDR), var=var(CDR))
}
# apply function to each combination of parameter values
library(plyr)
CDR_var <- mdply(parms, CDR_moments, zeta = zeta, 
                 iterations = iterations, stream_length = 600)
# plot results
library(ggplot2)
qplot(x = phi, y = var, color = gen_dist, linetype = gen_dist,
      data = CDR_var, geom = "point", alpha = I(0.5), size = I(1.5),
      xlab = "Prevalence", ylab = expression(Var(Y^C))) + 
  geom_smooth(method = "lm", formula = y ~ I(x^2 * (1 - x)^2) + 0, se = FALSE) + 
  labs(linetype="Generating distributions", color="Generating distributions") +
  theme_bw()
