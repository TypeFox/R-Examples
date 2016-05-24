
#' Mapping of Simulation Model Functions
#'
#' Provides information about the function for dynamic predictions associated
#' to a valid \code{simulation_model} key.
#' If \code{simulation_model} is missing or \code{NULL}, a character vector
#' of valid model keys is provided.
#' This function is designed as an assistant for using the functions
#' \code{\link{predict_inactivation}} and
#' \code{\link{fit_dynamic_inactivation}}.
#' For the adjustment of isothermal experiments with the function
#' \code{\link{fit_isothermal_inactivation}}, use the function
#' \code{\link{get_isothermal_model_data}}.
#'
#'
#' @param simulation_model (optional) character with a valid model key or
#'        \code{NULL}.
#'
#' @return If simulation_model is \code{NULL} or missing, a character vector of
#'      possible names. Otherwise, a list including information of the relevant
#'      function:
#'      \itemize{
#'          \item ode: Pointer to the function defining the model ode.
#'          \item cost: Pointer to the function calculating the error of the
#'                approximation.
#'          \item dtemp: logical defining whether the function requires the
#'                       definition of the first derivative of temperature.
#'          \item variables: a character vector defining which entry variables
#'                           are needed by the model.
#'          \item variables_priv: for internal use only.
#'          \item parameters: character vector with the parameters needed by
#'                            the model.
#'          }
#'
#' @export
#'
#' @seealso \code{\link{predict_inactivation}},
#'          \code{\link{fit_dynamic_inactivation}}
#'
get_model_data <- function(simulation_model = NULL) {

    function_map <- list(Bigelow = list(ode = dBigelow_model,
                                        dtemp = FALSE,
                                        variables = "N",
                                        variables_priv = "N",
                                        parameters = c("D_R", "z", "temp_ref")
                                        ),

#                          Bigelow_full = list(ode = dBigelow_model_full,
#                                                    dtemp = TRUE,
#                                                    variable = "N",
#                                                    parameters = c("D_R", "z",
#                                                                   "temp_ref")
#                                                    ),

                         Mafart = list(ode = dMafart_model,
                                             dtemp = FALSE,
                                             variables = "N",
                                             variables_priv = "N",
                                             parameters = c("delta_ref",
                                                            "temp_ref", "z",
                                                            "p")
                                             ),

#                          Weibull_Mafart_full = list(ode = dMafart_model_full,
#                                                     dtemp = TRUE,
#                                                     variable = "N",
#                                                     parameters = c("delta_ref",
#                                                                    "temp_ref",
#                                                                    "z", "p")
#                                                     ),

                         Geeraerd = list(ode = dGeeraerd_model,
                                         dtemp = FALSE,
                                         variables = c("N", "C_c"),
                                         variables_priv = c("N", "C_c"),
                                         parameters = c("D_R", "z", "N_min",
                                                        "temp_ref")
                                         ),

                         Peleg = list(ode = dPeleg_model,
                                            dtemp = FALSE,
                                            variables = "N",
                                            variables_priv = "logS",
                                            parameters = c("k_b", "temp_crit",
                                                           "n")
                                            )

#                          Weibull_Peleg_full = list(ode = dPeleg_model_full,
#                                                    dtemp = TRUE,
#                                                    variable = "N",
#                                                    parameters = c("k_b",
#                                                                   "temp_crit",
#                                                                   "n")
#                                                    )
                         )

    if (is.null(simulation_model)) {
        return(names(function_map))

    } else {

        function_data <- function_map[[simulation_model]]

        if (is.null(function_data)) {
            stop(paste("Unknown simulation model", simulation_model))

        } else {
            return(function_data)

        }
    }
}

#' Correctness Check of Model Parameters
#'
#' Makes 3 checks of compatibility between the input parameters for the
#' adjustment and the DOF of the inactivation model selected.
#' \itemize{
#'      \item Check of equal length of model DOF and input DOF. Raises a
#'            warning if failed.
#'      \item Check that every single one of the input DOF is a model DOF.
#'            Raises a warning if failed.
#'      \item Check that every single one of the model DOF are defined as
#'            an input DOF.}
#'
#' @param simulation_model character with a valid model key.
#' @param known_params named vector or list with the dof of the model
#'                     known.
#' @param starting_points named vector or list with the dof of the model to
#'                        be adjusted.
#' @param adjust_vars logical specifying whether the model variables need to
#'                    be included in the check (not used for isothermal
#'                    fit)
#'
check_model_params <- function(simulation_model, known_params,
                               starting_points, adjust_vars) {

    model_data <- get_model_data(simulation_model)

    variable_names <- paste0(model_data$variables, "0")

    if (adjust_vars) {
        model_dof <- c(model_data$parameters, variable_names)
    } else {
        model_dof <- c(model_data$parameters)
    }


    #- Check of lengths

    in_length <- length(known_params) + length(starting_points)

    if (length(model_dof) != in_length) {

        warning_msg <- paste0("Length of known_params + starting_points (",
                              in_length, ") not equal to model_dof (",
                              length(model_dof), ").")

        warning(warning_msg)
    }

    #- Check for input param names

    in_dof <- names(c(known_params, starting_points))

    for (each_dof in in_dof) {

        if (!(each_dof %in% model_dof)) {

            warning(paste0("Unknown input DoF: ", each_dof))
        }
    }

    #- Check for missing params

    for (each_dof in model_dof) {

        if (!(each_dof %in% in_dof)) {

            stop(paste("DOF ", each_dof, "not defined."))
        }
    }
}















