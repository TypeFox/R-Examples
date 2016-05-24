
#' Dynamic Prediction Intervals from a Monte Carlo Adjustment
#'
#' Given a model adjustment of a dynamic microbial inactivation process
#' performed using any of the functions in \code{bioinactivation} calculates
#' probability intervals at each time point using a Monte Carlo method.
#'
#' @param fit_object An object of classes \code{FitInactivationMCMC},
#'        \code{IsoFitInactivation} or \code{FitInactivation}.
#'
#' @param temp_profile data frame with discrete values of the temperature for
#'        each time. It must have one column named \code{time} and another named
#'        \code{temperature} providing discrete values of the temperature at
#'        time points.
#'
#' @param n_simulations a numeric indicating how many Monte Carlo simulations
#'        to perform. \code{100} by default.
#'
#' @param times numeric vector specifying the time points when results are
#'        desired. If \code{NULL}, the times in \code{MCMC_fit$best_prediction}
#'        are used. \code{NULL} by default.
#'
#' @param quantiles numeric vector indicating the quantiles to calculate in
#'        percentage. By default, it is set to c(2.5, 97.5) which generates a
#'        prediction interval with confidence 0.95. If \code{NULL}, the quantiles
#'        are not calculated and all the simulations are returned.
#'
#' @param additional_pars Additional parameters not included in the adjustment (e.g.
#'        the initial number of microorganism in an isothermal fit).
#'
#' @importFrom stats quantile median
#'
#' @return A data frame of class \code{PredInactivationMCMC}. On its first column,
#'         time at which the calculation has been made is indicated.
#'         If \code{quantiles = NULL}, the following columns contain the
#'         results of each simulation. Otherwise, the second and third columns
#'         provide the mean and median of the simulations at the given time
#'         point. Following columns contain the quantiles of the results.
#'
#' @export
#'
predict_inactivation_MCMC <- function(fit_object, temp_profile, n_simulations = 100,
                                      times = NULL, quantiles = c(2.5, 97.5),
                                      additional_pars = NULL){

    if (is.FitInactivationMCMC(fit_object)){

        MCMC_sample <- sample_MCMCfit(fit_object, times, n_simulations)

        par_sample <- MCMC_sample$par_sample
        simulation_model <- MCMC_sample$simulation_model
        known_pars = MCMC_sample$known_pars
        times <- MCMC_sample$times

    } else if (is.IsoFitInactivation(fit_object)){

        isoFit_sample <- sample_IsoFit(fit_object, times, n_simulations)

        par_sample <- isoFit_sample$par_sample
        simulation_model <- isoFit_sample$simulation_model
        known_pars = isoFit_sample$known_pars
        times <- isoFit_sample$times

    } else if (is.FitInactivation(fit_object)) {

        dynaFit_sample <- sample_dynaFit(fit_object, times, n_simulations)

        par_sample <- dynaFit_sample$par_sample
        simulation_model <- dynaFit_sample$simulation_model
        known_pars = dynaFit_sample$known_pars
        times <- dynaFit_sample$times

    } else {
        stop("fit_object must be a of class FitInactivationMCMC")
    }


    result_matrix <- matrix(nrow = length(times), ncol = n_simulations)

    for (i in 1:n_simulations){
        this_pars <- unlist(par_sample[i, ])
        this_prediction <- predict_inactivation(simulation_model,
                                                times,
                                                c(this_pars, known_pars, additional_pars),
                                                temp_profile
                                                )
        result_matrix[, i] <- this_prediction$simulation$N
    }

    ## Take means and quantiles from the results

    if (is.null(quantiles)) {
        out_matrix <- as.data.frame(cbind(times, result_matrix))

    } else {

        quantile_matrix <- apply(result_matrix, 1, quantile, probs = quantiles/100, na.rm = TRUE)
        median_matrix <- apply(result_matrix, 1, median, na.rm = TRUE)
        mean_matrix <- rowMeans(result_matrix, na.rm = TRUE)
        out_matrix <- as.data.frame(cbind(times,
                                          mean = mean_matrix,
                                          median = median_matrix,
                                          t(quantile_matrix)))

    }

    class(out_matrix) <- c("PredInactivationMCMC", class(out_matrix))
    return(out_matrix)
}

#'
#' Random sample of the parameters of a FitInactivationMCMC object
#'
#' Function to be called by \code{\link{predict_inactivation_MCMC}}. Produces
#' a random sample of the parameters calculated on the iterations of the Monte
#' Carlo simulation.
#'
#' @param MCMC_fit An object of class \code{FitInactivationMCMC} as generated
#'        by \code{\link{fit_inactivation_MCMC}}.
#' @param n_simulations a numeric indicating how many Monte Carlo simulations
#'        to perform.
#' @param times numeric vector specifying the time points when results are
#'        desired. If \code{NULL}, the times in \code{MCMC_fit$best_prediction}
#'        are used.
#'
#' @return Returns a list with 4 components.
#'         \itemize{
#'              \item par_sample: data frame with the parameter sample.
#'              \item simulation_model: model key for the simulation
#'              \item known_pars: parameters of the model known
#'              \item times: points where the calculation is made.
#'              }
#'
#' @importFrom dplyr sample_n
#'
sample_MCMCfit <- function(MCMC_fit, times, n_simulations){

    ## Extract fitted parameters and identify known pars

    pars <- as.data.frame(MCMC_fit$modMCMC$pars)
    fit_pars <- names(pars)

    best_pars <- MCMC_fit$best_prediction$model_parameters
    good_indexes <- !(names(best_pars) %in% fit_pars)
    known_pars <- unlist(best_pars[good_indexes])

    ## Extrat the rest of the model data

    simulation_model <- MCMC_fit$best_prediction$model

    if (is.null(times)) {
        times <- MCMC_fit$best_prediction$simulation$time
    }

#     temp_values <- MCMC_fit$best_prediction$temp_approximations$temp(times)
#     temp_profile <- data.frame(time = times,
#                                temperature = temp_values)

    ## Sample from pars and make the simulations

    par_sample <- sample_n(pars, n_simulations, replace = TRUE)

    ## Return the results

    out <- list(par_sample = par_sample,
                simulation_model = simulation_model,
                known_pars = known_pars,
                times = times
                )
    return(out)
}

#'
#' Random sample of the parameters of a IsoFitInactivation object
#'
#' Function to be called by \code{\link{predict_inactivation_MCMC}}. Produces
#' a random sample of the parameters adjusted from isothermal experiments.
#'
#' It is assumed that the parameters follow a Multivariate Normal distribution
#' with the mean and covariance matrix estimated by the adjustment. The function
#' produces a random sample using the function \code{\link{mvrnorm}}.
#'
#' @param iso_fit An object of class \code{FitInactivationMCMC} as generated
#'        by \code{\link{fit_inactivation_MCMC}}.
#' @param n_simulations a numeric indicating how many Monte Carlo simulations
#'        to perform.
#' @param times numeric vector specifying the time points when results are
#'        desired. If \code{NULL}, an equispaced interval between \code{0}
#'        and the maximum time of the observations with length \code{50}
#'        is used.
#'
#' @return Returns a list with 4 components.
#'         \itemize{
#'              \item par_sample: data frame with the parameter sample.
#'              \item simulation_model: model key for the simulation
#'              \item known_pars: parameters of the model known
#'              \item times: points where the calculation is made.
#'              }
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats coef
#' @importFrom stats vcov
#'
#' @seealso \code{\link{mvrnorm}}
#'
sample_IsoFit <- function(iso_fit, times, n_simulations){

    ## Take the sample of the parameters

    par_sample <- as.data.frame(mvrnorm(n_simulations,
                                        coef(iso_fit$nls),
                                        vcov(iso_fit$nls))
                                )

    ## Remove rows where any parameter is below 0

    bad_index <- apply(par_sample, 1, function(x) any(x<0))
    par_sample <- par_sample[!bad_index, ]

    ## Identify known pars

    fit_pars <- names(par_sample)

    best_pars <- iso_fit$parameters
    good_indexes <- !(names(best_pars) %in% fit_pars)
    known_pars <- unlist(best_pars[good_indexes])

    ## Extrat the rest of the model data

    simulation_model <- iso_fit$model

    if (is.null(times)) {
        times <- seq(0, max(iso_fit$data$time), length = 50)
    }

    ## Return the results

    out <- list(par_sample = as.data.frame(par_sample),
                simulation_model = simulation_model,
                known_pars = known_pars,
                times = times
    )
    return(out)

}

#'
#' Random sample of the parameters of a FitInactivation object
#'
#' Function to be called by \code{\link{predict_inactivation_MCMC}}. Produces
#' a random sample of the parameters adjusted from dynamic experiments with non
#' linear regression.
#'
#' It is assumed that the parameters follow a Multivariate Normal distribution
#' with the mean the parameters estimated by \code{\link{modFit}}. The unscaled
#' covariance matrix returned by \code{\link{modFit}} is used.
#'
#' The function produces a random sample using the function \code{\link{mvrnorm}}.
#'
#' @param dynamic_fit An object of class \code{FitInactivationMCMC} as generated
#'        by \code{\link{fit_inactivation_MCMC}}.
#' @param n_simulations a numeric indicating how many Monte Carlo simulations
#'        to perform.
#' @param times numeric vector specifying the time points when results are
#'        desired. If \code{NULL}, the times in \code{MCMC_fit$best_prediction}
#'        are used.
#'
#' @return Returns a list with 4 components.
#'         \itemize{
#'              \item par_sample: data frame with the parameter sample.
#'              \item simulation_model: model key for the simulation
#'              \item known_pars: parameters of the model known
#'              \item times: points where the calculation is made.
#'              }
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats coef
#' @importFrom stats vcov
#'
sample_dynaFit <- function(dynamic_fit, times, n_simulations){

    ## Take the sample of the parameters

    sP <- summary(dynamic_fit)

    par_sample <- as.data.frame(mvrnorm(n_simulations,
                                        coef(dynamic_fit$fit_results),
                                        sP$cov.unscaled)
                                )

    ## Remove rows where any parameter is below 0

    bad_index <- apply(par_sample, 1, function(x) any(x<0))
    par_sample <- par_sample[!bad_index, ]

    ## Identify known pars

    fit_pars <- names(par_sample)

    best_pars <- dynamic_fit$best_prediction$model_parameters
    good_indexes <- !(names(best_pars) %in% fit_pars)
    known_pars <- unlist(best_pars[good_indexes])

    ## Extrat the rest of the model data

    simulation_model <- dynamic_fit$best_prediction$model

    if (is.null(times)) {
        times <- seq(0, max(dynamic_fit$data$time), length = 50)
    }

    ## Return the results

    out <- list(par_sample = as.data.frame(par_sample),
                simulation_model = simulation_model,
                known_pars = known_pars,
                times = times
                )
    return(out)

}






