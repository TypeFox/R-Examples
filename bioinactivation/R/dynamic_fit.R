
#' Fitting of Dynamic Inactivation Models
#'
#' Fits the parameters of an inactivation model to experimental data.
#' The function \code{\link{modFit}} of the package \code{\link{FME}} is
#' used for the adjustment.
#'
#' @param experiment_data data frame with the experimental data to be adjusted.
#'        It must have a column named \dQuote{time} and another one named
#'        \dQuote{N}.
#' @param simulation_model character identifying the model to be used.
#' @param temp_profile data frame with discrete values of the temperature for
#'        each time. It must have one column named \code{time} and another named
#'        \code{temperature} providing discrete values of the temperature at
#'        time points.
#' @param starting_points starting values of the parameters for the adjustment.
#' @param upper_bounds named numerical vector defining the upper bounds of the
#'        parameters for the adjustment.
#' @param lower_bounds named numerical vector defining the lower bounds of the
#'        parameters for the adjustment.
#' @param known_params named numerical vector with the fixed (i.e., not
#'        adjustable) model parameters.
#' @param minimize_log logical. If \code{TRUE}, the adjustment is based on the
#'        minimization of the error of the logarithmic count. \code{TRUE} by
#'        default.
#' @param tol0 numeric. Observations at time 0 make Weibull-based models singular.
#'        The time for observatins taken at time 0 are changed for this value.
#' @param ... further arguments passed to \code{\link{modFit}}
#'
#' @importFrom dplyr mutate_
#' @importFrom dplyr %>%
#' @importFrom dplyr select_
#' @importFrom FME modFit
#' @importFrom lazyeval interp
#'
#' @return A list of class \code{FitInactivation} with the following items:
#'      \itemize{
#'          \item fit_results: a list of class \code{modFit} with the info
#'                of the adjustment.
#'          \item best_prediction: a list of class \code{SimulInactivation}
#'                with prediction made by the adjusted model.
#'          \item data: a data frame with the data used for the fitting.
#'          }
#'
#' @export
#'
#' @seealso \code{\link{modFit}}
#'
#' @examples
#' ## EXAMPLE 1 ------
#'
#' data(dynamic_inactivation)  # The example data set is used.
#'
#' get_model_data()  # Retrieve the valid model keys.
#'
#' simulation_model <- "Peleg"  # Peleg's model will be used
#'
#' model_data <- get_model_data(simulation_model)
#' model_data$parameters  # Set the model parameters
#'
#' dummy_temp <- data.frame(time = c(0, 1.25, 2.25, 4.6),
#'                          temperature = c(70, 105, 105, 70))  # Dummy temp. profile
#'
#' ## Set known parameters and initial points/bounds for unknown ones
#'
#' known_params = c(temp_crit = 100)
#'
#' starting_points <- c(n = 1, k_b = 0.25, N0 = 1e+05)
#' upper_bounds <- c(n = 2, k_b = 1, N0 = Inf)
#' lower_bounds <- c(n = 0, k_b = 0, N0 = 1e4)
#'
#' dynamic_fit <- fit_dynamic_inactivation(dynamic_inactivation, simulation_model,
#'                                         dummy_temp, starting_points,
#'                                         upper_bounds, lower_bounds,
#'                                         known_params)
#'
#' plot(dynamic_fit)
#'
#' ## END EXAMPLE 1 -----
#'
fit_dynamic_inactivation <- function(experiment_data, simulation_model, temp_profile,
                                     starting_points, upper_bounds, lower_bounds,
                                     known_params, ..., minimize_log = TRUE, tol0 = 1e-5) {

    #- Check of the model parameters

    check_model_params(simulation_model, known_params, starting_points, TRUE)

    #- Gather the information

    model_data <- get_model_data(simulation_model)

    if (minimize_log) {

        data_for_fit <- mutate_(experiment_data,
                                logN = interp(~log10(N), N=as.name("N")))

        data_for_fit <- select_(data_for_fit, as.name("time"), as.name("logN"))

    } else {
        data_for_fit <- select_(experiment_data, as.name("time"), as.name("N"))
    }

    #- Add a small tolerance to data at time 0 to avoid singularities

    data_0 <- data_for_fit$time < tol0
    data_for_fit[data_0, "time"] <- tol0

    #- Call the fitting function

    Fit <- modFit(f = get_prediction_cost, p = starting_points,
                  temp_profile = temp_profile,
                  data_for_fit = data_for_fit,
                  lower = lower_bounds, upper = upper_bounds,
                  known_params = known_params,
                  simulation_model = simulation_model,
                  ...
                  )

    #- Output the results

    pars_fit <- Fit$par

    prediction_time <- seq(min(data_for_fit$time), max(data_for_fit$time), length=100)
    best_prediction <- predict_inactivation(simulation_model,
                                            prediction_time,
                                            as.list(c(pars_fit, known_params)),
                                            temp_profile
                                            )

    out_results <- list(fit_results = Fit,
                        best_prediction = best_prediction,
                        data = data_for_fit
                        )

    class(out_results) <- c("FitInactivation", class(out_results))
    return(out_results)

}






