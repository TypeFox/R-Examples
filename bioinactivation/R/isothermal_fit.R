
#' Fit of Isothermal Experiments
#'
#' Fits the parameters of the model chosen to a set of isothermal experiments
#' using nonlinear regression through the function \code{\link{nls}}.
#'
#' @param model_name character specyfing the model to adjust.
#' @param death_data data frame with the experiment data where each row is one
#'        observation. It must have the following columns:
#'        \itemize{
#'          \item log_diff: Number of logarithmic reductions at each data point.
#'          \item temp: Temperature of the data point.
#'          \item time: Time of the data point.
#'          }
#' @param starting_point List with the initial values of the parameters
#'        for the adjustment.
#' @param known_params List of the parameters of the model known.
#' @param adjust_log logical. If \code{TRUE}, the adjustment is based on the
#'        minimization of the error of the logarithmic microbial count. If
#'        \code{FALSE}, it is based on the minimization of the error of the
#'        microbial count. \code{TRUE} by default.
#'
#' @return An instance of class \code{IsoFitInactivation} with the results.
#'         This list has four entries:
#'         \itemize{
#'           \item nls: The object of class \code{\link{nls}} with the results
#'                      of the adjustment.
#'           \item parameters: a list with the values of the model parameters,
#'                             both fixed and adjusted.
#'           \item model: a string with the key identifying the model used.
#'           \item data: the inactivation data used for the fit.
#'           }
#'
#' @export
#'
#' @importFrom stats coefficients
#'
#' @seealso \code{\link{nls}}
#'
#' @examples
#' ## EXAMPLE 1 -----------
#'
#' data(isothermal_inactivation)  # data set used for the example.
#'
#' get_isothermal_model_data()  # retrieve valid model keys.
#' model_name <- "Bigelow"  # Bigelow's model will be used for the adjustment.
#'
#' model_data <- get_isothermal_model_data(model_name)
#' model_data$params  # Get the parameters of the model
#'
#' ## Define the input arguments
#'
#' known_params = list(temp_ref = 100)
#' starting_point <- c(z = 10,D_R = 1)
#' adjust_log <- TRUE
#'
#' ## Call the fitting function
#' iso_fit <- fit_isothermal_inactivation(model_name,
#'                                        isothermal_inactivation, starting_point,
#'                                        known_params)
#'
#' ## Output of the results
#'
#' plot(iso_fit)
#'
#' ## END EXAMPLE 1 --------
#'
fit_isothermal_inactivation <- function(model_name, death_data, starting_point,
                                        known_params, adjust_log = TRUE) {

    #- Check of the model parameters

    check_model_params(model_name, known_params, starting_point, FALSE)

    #- Make the adjustment

    adjust_results <- with(known_params, {

        model_data <- get_isothermal_model_data(model_name)

        if (is.null(model_data)) {
            stop(paste("Unknown model:", model_data))
        }

        if (adjust_log) {

            my_formula = as.formula( paste("log_diff ~", model_data$formula_iso) )

        } else {
            death_data <- mutate(death_data, S = 10^log_diff)
            my_formula = as.formula( paste("S ~10^", model_data$formula_iso) )

        }

        adjust_results <- nls(my_formula,
                              data = death_data,
                              start = starting_point)

        return(adjust_results)

    })

    #- Output of the results

    out <- list()
    class(out) <- c("IsoFitInactivation", class(out))

    out$nls <- adjust_results
    out$parameters <- c(known_params, coefficients(adjust_results))
    out$model <- model_name
    out$data <- death_data

    return(out)

}

#' Isothermal Model Data
#'
#' Provides information of the models implemented for fitting of isothermal
#' data.
#' This models are valid only for isothermal adjustment with the function
#' \code{\link{fit_isothermal_inactivation}}. To make predictions with the
#' function \code{\link{predict_inactivation}} or adjust dynamic experiments
#' with \code{\link{fit_dynamic_inactivation}}, use
#' \code{\link{get_model_data}}.
#'
#' @param model_name Optional string with the key of the model to use.
#' @return If \code{model_name} is missing, a list of the valid model keys.
#'         If \code{model_name} is not a valid key, NULL is returned.
#'         Otherwise, a list with the parameters of the model selected and its
#'         \code{formula} for the nonlinear adjustment.
#'
#' @export
#'
get_isothermal_model_data <- function(model_name = "valids") {

    switch(as.character(model_name),
           Mafart = list(params = c("delta_ref", "z", "p", "temp_ref"),
                         formula_iso = "WeibullMafart_iso(time, temp, delta_ref, z, p, temp_ref)",
                         prediction = "WeibullMafart_iso"
                         ),

           Peleg = list(params = c("n", "k_b", "temp_crit"),
                        formula_iso = "WeibullPeleg_iso(time, temp, n, k_b, temp_crit)",
                        prediction = "WeibullPeleg_iso"
                        ),

           Bigelow = list(params = c("z", "D_R", "temp_ref"),
                        formula_iso = "Bigelow_iso(time, temp, z, D_R, temp_ref)",
                        prediction = "Bigelow_iso"
                        ),

           valids = c("Mafart", "Peleg", "Bigelow")
           )

}
