x <- matrix(runif(50), 10, 5)
beta <- matrix(1, 10, 5)
gamma <- matrix(.5, 10, 5)

ret <- batch_normalization(x, gamma, beta)

mu <- ret[[1]]
sigma_2 <- ret[[2]]
x_hat <- ret[[3]]
y <- ret[[4]]

target <- matrix(1, 10, 5)
delta_y <- y - target

ret <- batch_normalization_differential(delta_y,
                                        mu,
                                        sigma_2,
                                        x,
                                        x_hat,
                                        y,
                                        gamma,
                                        beta)

delta_x <- ret[[1]]
delta_gamma <- ret[[2]]
delta_beta <- ret[[3]]
delta_x_hat <- ret[[4]]
delta_sigma_2 <- ret[[5]]
delta_mu <- ret[[6]]

write_2_csv <- function(data, file_name) {
  file_name <- paste0('excl/test_batch_normalization_differential/', file_name, '.csv')
  write.csv(data, file = file_name)
}

write_2_csv(x, "x.csv")
write_2_csv(y, "y.csv")
write_2_csv(delta_y, "delta_y.csv")
write_2_csv(delta_x, "delta_x.csv")
write_2_csv(delta_gamma, "delta_gamma.csv")
write_2_csv(delta_beta, "detla_beta.csv")

write_2_csv(delta_x_hat, "delta_x_hat")
write_2_csv(delta_sigma_2, "detla_sigma_2")
write_2_csv(delta_mu, "delta_mu")

write_2_csv(mu, "mu")
write_2_csv(sigma_2, "sigma_2")

