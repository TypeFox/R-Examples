# # Unit tests for bayes weights
#
# #source("./R/bayes_weights.R")
#
# # Plots  ----------------
# # generate means
# J <- 2000
# mu <- -abs(rnorm(J))
# sigma <- 1 * rep(1, J)
# alpha <- 1 / J
#
# # find weights
# w_1 <- bayes_weights(mu, sigma, alpha)
#
# plot(mu, w_1$w)
#
# # vary alpha
# for (j in (1:4)) {
#     alpha <- 2^j / J
#     w_1 <- bayes_weights(mu, sigma, alpha)
#     plot(mu, w_1$w)
# }
#
# # Brent  ----------------
# # generate means
# J <- 200
# mu <- -abs(rnorm(J))
# sigma <- 1 * rep(1, J)
# alpha <- 100
#
# # find weights
# #source("bayes_weights.R")
# w_1 <- bayes_weights(mu, sigma, q = alpha/J)
# plot(mu,w_1$w)
#
