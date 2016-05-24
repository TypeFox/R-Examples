
#' Error of the Prediction of Microbial Inactivation
#'
#' Calculates the error of the prediction of microbial inactivation for
#' the chosen inactivation model and the given parameters with respect to
#' the experimental data provided.
#' This function is compatible with the function
#' \code{\link{fit_dynamic_inactivation}}.
#'
#' @param data_for_fit A data frame with the experimental data to fit. It
#'        must contain a column named \dQuote{time} and another one named
#'        \dQuote{N}.
#' @param temp_profile \code{data.frame} defining the temperature profile. It
#'        must have a column named \dQuote{time} and another named
#'        \dQuote{temperature}.
#' @param simulation_model character key defining the inactivation model.
#' @param P list with the unknown parameters of the model to be adjusted.
#' @param known_params list with the parameters of the model fixed (i.e.,
#'        not adjusted)
#'
#' @importFrom FME modCost
#' @importFrom stats complete.cases
#'
#' @return An instance of \code{\link{modCost}} with the error of the
#'         prediction.
#'
#' @seealso \code{\link{modCost}}, \code{\link{fit_dynamic_inactivation}}
#'
get_prediction_cost <- function(data_for_fit, temp_profile,
                                simulation_model,
                                P, known_params
                                ) {

    prediction_data <- predict_inactivation(simulation_model = simulation_model,
                                            times = sort(unique(data_for_fit$time)),
                                            parms = as.list(c(P, known_params)),
                                            temp_profile = temp_profile
                                            )

    prediction <- prediction_data$simulation
    prediction <- prediction[names(data_for_fit)]  # Take only the relevant columns
    prediction <- prediction[complete.cases(prediction),]  # NAs produced by negative values in x. This is caused by numerical error when N is very small.

    cost <- modCost(model = prediction, obs = data_for_fit)
    return(cost)

}
