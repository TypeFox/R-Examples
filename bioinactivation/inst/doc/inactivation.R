## ------------------------------------------------------------------------
library(bioinactivation)

## ------------------------------------------------------------------------
data(isothermal_inactivation)

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
library(ggplot2)
ggplot(data = isothermal_inactivation) +
    geom_point(aes(x = time, y = log_diff, col = as.factor(temp))) +
    ggtitle("Example dataset: isothermal_inactivation")


## ------------------------------------------------------------------------
data(dynamic_inactivation)

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
ggplot(data = dynamic_inactivation) +
    geom_point(aes(x = time, y = log10(N))) +
    ggtitle("Example dataset: dynamic_inactivation. Observations")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
ggplot(data = dynamic_inactivation) +
    geom_point(aes(x = time, y = temperature))  +
    ggtitle("Example dataset: dynamic_inactivation. Temperature profile")

## ------------------------------------------------------------------------
data(laterosporus_iso)

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
ggplot(data = laterosporus_iso) +
    geom_point(aes(x = time, y = log_diff, col = as.factor(temp))) +
    ggtitle("Example dataset: laterosporus_iso")


## ------------------------------------------------------------------------
data(laterosporus_dyna)

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
ggplot(data = laterosporus_dyna) +
    geom_point(aes(x = time, y = logN)) +
    ggtitle("Example dataset: laterosporus_dyna. Observations")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
ggplot(data = laterosporus_dyna) +
    geom_point(aes(x = time, y = temp))  +
    ggtitle("Example dataset: laterosporus_dyna. Temperature profile")

## ----, echo=FALSE, fig.width=6, fig.height=4, fig.align='center'---------
Temperature <- seq(30, 80, length=100)
b <- log( 1 + exp(0.3*(Temperature - 60)))
b_1 <- log( 1 + exp(0.2*(Temperature - 60)))
plot(Temperature, b, type = "l", col="red")
lines(Temperature, b_1, type = "l", col = "blue")
legend("topleft", c("k = 0.3","k = 0.2"), col=c("red", "blue"), lwd=1, cex = 1)

## ------------------------------------------------------------------------
get_model_data()

## ------------------------------------------------------------------------
example_model <- "Geeraerd"

## ------------------------------------------------------------------------
times <- seq(0, 5, length=100)

## ------------------------------------------------------------------------
model_data <- get_model_data(example_model)
print(model_data$parameters)
print(model_data$variables)

## ------------------------------------------------------------------------
model_parms <- c(D_R = 1,
                 z = 10,
                 N_min = 1e2,
                 temp_ref = 100,
                 N0 = 1e5,
                 C_c0 = 1e1
                 )

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
temperature_profile <- data.frame(time = c(0, 5),
                                  temperature = c(70, 120))
plot(temperature_profile, type = "l")
title("Example temperature profile")

## ------------------------------------------------------------------------
prediction_results <- predict_inactivation(example_model, times, model_parms, temperature_profile)

## ------------------------------------------------------------------------
head(prediction_results$simulation)

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
p <- plot(prediction_results)
print(p)

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
library(ggplot2)
p + theme_light() + xlab("time (min)")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
plot(prediction_results, make_gg = FALSE,
     xlab = "Time (min)", ylab = "logN")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
parms_no_shoulder <- c(D_R = 1,
                       z = 10,
                       N_min = 100,
                       temp_ref = 100,
                       N0 = 100000,
                       C_c0 = 0
                       )


prediction_no_shoulder <- predict_inactivation(example_model, times, parms_no_shoulder, temperature_profile)

plot(prediction_no_shoulder)


## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
parms_no_tail <- c(D_R = 1,
                   z = 10,
                   N_min = 0,
                   temp_ref = 100,
                   N0 = 100000,
                   C_c0 = 100
                   )


prediction_no_tail <- predict_inactivation(example_model, times, parms_no_tail, temperature_profile)

plot(prediction_no_tail)


## ------------------------------------------------------------------------
data(isothermal_inactivation)

## ------------------------------------------------------------------------
get_isothermal_model_data()

## ------------------------------------------------------------------------
model_name <- "Bigelow"

## ------------------------------------------------------------------------
model_data <- get_isothermal_model_data(model_name)
model_data$params

## ------------------------------------------------------------------------
known_params = list(temp_ref = 100)

## ------------------------------------------------------------------------
starting_point <- list(z = 10,
                       D_R = 1.5
                       )

## ------------------------------------------------------------------------
iso_fit <- fit_isothermal_inactivation(model_name, isothermal_inactivation,  
                                       starting_point, known_params)

## ------------------------------------------------------------------------
summary(iso_fit$nls)
vcov(iso_fit$nls)  # Calculates variance-covariance matrix {stats}
confint(iso_fit$nls)  # Calculates confidence intervals {stats}

## ------------------------------------------------------------------------
iso_fit$parameters

## ------------------------------------------------------------------------
iso_fit$model

## ------------------------------------------------------------------------
head(iso_fit$data)

## ------------------------------------------------------------------------
plot(iso_fit, ylab = "Number of logarithmic reductions",
     xlab = "Time (min)")

## ------------------------------------------------------------------------
data(dynamic_inactivation)

## ------------------------------------------------------------------------
get_model_data()

## ------------------------------------------------------------------------
simulation_model <- "Peleg"

## ------------------------------------------------------------------------
dummy_temp <- data.frame(time = dynamic_inactivation$time,
                         temperature = dynamic_inactivation$temperature)

## ------------------------------------------------------------------------
model_data <- get_model_data(simulation_model)
model_data$parameters
model_data$variables

## ------------------------------------------------------------------------
known_params = c(temp_crit = 100)

## ------------------------------------------------------------------------
starting_points <- c(n = 1,
                     k_b = 0.25,
                     N0 = 1e+05)
upper_bounds <- c(n = 2,
                  k_b = 1,
                  N0 = Inf)

lower_bounds <- c(n = 0,
                  k_b = 0,
                  N0 = 1e4)

## ------------------------------------------------------------------------
dynamic_fit <- fit_dynamic_inactivation(dynamic_inactivation, simulation_model, dummy_temp,  
                                        starting_points, upper_bounds, lower_bounds,  
                                        known_params)

## ------------------------------------------------------------------------
summary(dynamic_fit$fit_results)

## ------------------------------------------------------------------------
head(dynamic_fit$data)

## ------------------------------------------------------------------------
dynamic_fit$best_prediction$model_parameters

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
p <- plot(dynamic_fit)
p + theme_light()

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
plot(dynamic_fit, make_gg = FALSE,
     xlab = "Time (min)", ylab = "logN")

## ------------------------------------------------------------------------
set.seed(82619)

MCMC_fit <- fit_inactivation_MCMC(dynamic_inactivation, simulation_model, dummy_temp,  
                                        starting_points, upper_bounds, lower_bounds,  
                                        known_params, niter = 50)

## ------------------------------------------------------------------------
summary(MCMC_fit$modMCMC)

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
p <- plot(MCMC_fit)
p + xlab("Time (min)")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
plot(MCMC_fit, make_gg = FALSE,
     xlab = "Time (min)", ylab = "logN")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
library(MASS)
sP <- summary(dynamic_fit)
sample_iso <- mvrnorm(1000, coef(dynamic_fit$fit_results), sP$cov.unscaled)

ggplot(as.data.frame(sample_iso)) +
    geom_point(aes(x = N0, y = k_b), alpha = 0.5, size = 3) +
    ggtitle("Example of multivariate normal sampling")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
dummy_temp <- data.frame(time = c(0, 1, 4, 6), temperature = c(80, 110, 110, 80))
ggplot(dummy_temp) +
    geom_line(aes(x = time, y = temperature)) +
    ggtitle("Temperature profile used for the prediction interval")

## ------------------------------------------------------------------------
set.seed(7163)
pred_MCMC_1 <- predict_inactivation_MCMC(dynamic_fit, dummy_temp)
head(pred_MCMC_1)

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
p <- plot(pred_MCMC_1)
p + xlab("time")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
set.seed(7163)
pred_MCMC_2 <- predict_inactivation_MCMC(dynamic_fit, dummy_temp, quantiles = NULL)
plot(pred_MCMC_2, make_gg = FALSE,
     xlab = "Time (min)", ylab = "logN")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
set.seed(7163)
pred_MCMC_3 <- predict_inactivation_MCMC(dynamic_fit, dummy_temp, n_simulations = 200)
p <- plot(pred_MCMC_3)
p + ylab("logN")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
set.seed(7163)
pred_MCMC_4 <- predict_inactivation_MCMC(dynamic_fit, dummy_temp, quantiles = c(0, 95))
p <- plot(pred_MCMC_4)
p + ylab("logN")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
library(MASS)
sample_iso <- mvrnorm(1000, coef(iso_fit$nls), vcov(iso_fit$nls))

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
set.seed(918734)
times <- seq(0, 6, length = 50)
dummy_temp <- data.frame(time = c(0, 3, 4, 6), temperature = c(80, 110, 110, 80))
pred_int_iso <- predict_inactivation_MCMC(iso_fit, dummy_temp,
                                          times = times,
                                          additional_pars = c(N0 = 1e5))
p <- plot(pred_int_iso)
p + xlab("time")

## ----, fig.width=6, fig.height=4, fig.align='center'---------------------
set.seed(918734)
dummy_temp <- data.frame(time = c(0, 1, 4, 6), temperature = c(80, 110, 110, 80))
pred_int_MCMC <- predict_inactivation_MCMC(MCMC_fit, dummy_temp)
p <- plot(pred_int_MCMC)
p + xlab("time")

