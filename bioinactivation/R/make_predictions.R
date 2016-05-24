
#' Prediction of Dynamic Inactivation
#'
#' Predicts the inactivation of a microorganism under isothermal or
#' non-isothermal temperature conditions. The thermal resistence of
#' the microorganism are defined with the input arguments.
#'
#' The value of the temperature is calculated at each value of time by
#' linear interpolation of the values provided by the input argument
#' \code{temp_profile}.
#' The function \code{\link{ode}} of the package \code{\link{deSolve}} is
#' used for the resolution of the differential equation.
#'
#' @param simulation_model character identifying the model to be used.
#' @param times numeric vector of output times.
#' @param parms list of parameters defining the parameters of the model.
#' @param temp_profile data frame with discrete values of the temperature for
#'        each time. It must have one column named \code{time} and another named
#'        \code{temperature} providing discrete values of the temperature at
#'        time points.
#' @param tol0 numeric. Observations at time 0 make Weibull-based models singular.
#'        The time for observatins taken at time 0 are changed for this value.
#'        By default (`tol0 = 1e-5`)
#'
#' @importFrom deSolve ode
#' @importFrom dplyr mutate_
#' @importFrom lazyeval interp
#'
#' @return A list of class \code{SimulInactivation} with the results. It has
#'         the following entries:
#'         \itemize{
#'           \item model: character defining the model use for the prediction.
#'           \item model_parameters: named numeric vector with the values of
#'                                   the model parameters used.
#'           \item temp_approximations: function used for the interpolation of
#'                                      the temperature. For a numeric value of
#'                                      time given, returns the value of the
#'                                      temperature and its first derivative.
#'           \item simulation: A data frame with the results calculated. Its
#'                             first column contains the times at which the
#'                             solution has been calculated. The following
#'                             columns the values of the variables of the
#'                             model. The three last columns provide the
#'                             values of \code{logN}, \code{S} and
#'                             \code{logS}.
#'           }
#'
#'
#' @export
#'
#' @seealso \code{\link{ode}}, \code{\link{get_model_data}}
#'
#' @examples
#' ## EXAMPLE 1 -----------
#'
#' ## Retrieve the model keys available for dynamic models.
#' get_model_data()
#'
#' ## Set the input arguments
#' example_model <- "Geeraerd"  # Geeraerd's model will be used
#' times <- seq(0, 5, length=100)  # values of time for output
#'
#' model_data <- get_model_data(example_model)  # Retrive the data of the model used
#' print(model_data$parameters)
#' print(model_data$variables)
#' model_parms <- c(D_R = 1,
#'                  z = 10,
#'                  N_min = 100,
#'                  temp_ref = 100,
#'                  N0 = 100000,
#'                  C_c0 = 1000
#'                  )
#'
#' ## Define the temperature profile for the prediction
#' temperature_profile <- data.frame(time = c(0, 1.25, 2.25, 4.6),
#'                                   temperature = c(70, 105, 105, 70))
#'
#' ## Call the prediction function
#' prediction_results <- predict_inactivation(example_model, times,
#'                                            model_parms, temperature_profile)
#'
#' ## Show the results
#' head(prediction_results$simulation)
#' plot(prediction_results)
#'
#' ## END EXAMPLE 1 -----------
#'
predict_inactivation <- function(simulation_model, times, parms, temp_profile, tol0 = 1e-5){

    #- Check of the model parameters

    check_model_params(simulation_model, c(), parms, TRUE)

    #- Remove 0 values of times

    times[times < tol0] <- tol0

    #- Gather the information

    model_data <- get_model_data(simulation_model)
    temp_approximations <- build_temperature_interpolator(temp_profile)

    #- Build the initial value

    if (all(model_data$variables_priv == "N")) {

        xstart <- c(N = parms[["N0"]])

    } else if (all(model_data$variables_priv == c("N", "C_c"))) {


        xstart <- c(N = parms[["N0"]], C_c = parms[["C_c0"]])

    } else if (all(model_data$variables_priv == "logS")) {

        xstart <- c(logS = -1e-6)
    }

    #- Call the ode function

    if (model_data$dtemp) {

        out <- ode(y = xstart, times = sort(times),
                   func = model_data$ode, parms = parms,
                   temp_profile = temp_approximations$temp,
                   dtemp_profile = temp_approximations$dtemp)

    } else {

        out <- ode(y = xstart, times = sort(times),
                   func = model_data$ode, parms = parms,
                   temp_profile = temp_approximations$temp)

    }

    #- Prepare result for output

    out_value <- list()
    out_value$model <- simulation_model
    out_value$model_parameters <- parms
    out_value$temp_approximations <- temp_approximations
    out <- as.data.frame(out)

    if ("N" %in% names(out)) {

        out <- mutate_(out,
                       logN = interp(~log10(N), N=as.name("N")),
                       S = interp(~N/parms[["N0"]], N=as.name("N")))

        out <- mutate_(out, logS = interp(~log10(S), S=as.name("S")))

    } else {

        out <- mutate_(out,
                       S = interp(~10^logS, logS=as.name("logS")),
                       N = interp(~S*parms[["N0"]], S=as.name("S")))

        out <- mutate_(out, logN = interp(~log10(N), N=as.name("N")))

    }

    out_value$simulation <- out

    class(out_value) <- c("SimulInactivation", class(out_value))

    return(out_value)

}


