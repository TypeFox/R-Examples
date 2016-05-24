
#'
#' First Derivate of the Linear Bigelow Model
#'
#' Calculates the first derivative of the linearized version of Bigelow's model
#' for dynamic problems at a given time for the model parameters provided and
#' the environmental conditions given.
#'
#' The model is developed from the isothermal Bigelow's model assuming during
#' the derivation process that \eqn{D_T} is time independenent.
#'
#' This function is compatible with the function
#' \code{\link{predict_inactivation}}.
#'
#' @section Model Equation:
#'
#'    \deqn{\frac{dN}{dt} = - N \frac {\mathrm{ln}(10)}{D_T(T)}}{
#'          dN/dt = - N * ln(10)/D_T(T)}
#'
#' @param t numeric vector indicating the time of the experiment.
#' @param x list with the value of N at t.
#' @param parms parameters for the secondary model. No explicit check of their
#'        validity is performed.
#' @param temp_profile a function that provides the temperature at a given
#'        time.
#'
#' @return The value of the first derivative of \eqn{N} at time \code{t} as a
#'         list.
#'
#' @section Model parameters:
#'      \itemize{
#'          \item temp_ref: Reference temperature for the calculation,
#'          \item D_R: D-value at the reference temperature,
#'          \item z: z value.
#'          }
#'
#' @seealso \code{\link{predict_inactivation}}
#'
dBigelow_model <- function(t, x, parms, temp_profile)  {

    with(as.list(c(parms, x)), {

        temp <- temp_profile(t)

        D_T <- D_R * 10^( (temp_ref - temp)/z )

        dN <- - N * log(10)/D_T

        res <- c(dN)
        return(list(res))

    })
}


# #'
# #' First Derivate of the Bigelow's Model with Full Derivation
# #'
# #' Calculates the first derivative of Bigelow's model at a given time for
# #' the model parameters provided and the environmental conditions given.
# #'
# #' The model is developed from the isothermal Bigelow's model taking into
# #' account in the derivation the time dependence of \eqn{D_T} for
# #' non-isothermal temperature profiles.
# #' For the linearized version, see \code{\link{dBigelow_model_linear}}
# #'
# #' This function is compatible with the function
# #' \code{\link{predict_inactivation}}.
# #'
# #' @section Model Equation:
# #'
# #' \deqn{\frac{dN}{dt} = - N \frac{ \mathrm{ln}(10) }{D_T}
# #'                      (1 + \mathrm{ln}(10) \frac{t}{z} \frac{dT}{dt}}{
# #'     dN/dt = - N * ln(10) / D_T(T) * (1 + ln(10)*t/z * dT/dt)}
# #'
# #' @param t numeric vector indicating the time of the experiment.
# #' @param x list with the value of N at t.
# #' @param parms parameters for the secondary model. No explicit check of their
# #'        validity is performed (see section \bold{Model Parameters}).
# #' @param temp_profile a function that provides the temperature at a given time.
# #' @param dtemp_profile a function that provides the derivative of the
# #'        temperature at a given time.
# #'
# #' @return The value of the first derivative of \eqn{N} at time \code{t} as a
# #'         list.
# #'
# #' @section Model parameters:
# #'      \itemize{
# #'          \item temp_ref: Reference temperature for the calculation,
# #'          \item D_R: D-value at the reference temperature,
# #'          \item z: z value.
# #'          }
# #'
# #' @seealso \code{\link{predict_inactivation}}
# #'
# dBigelow_model_full <- function(t, x, parms, temp_profile, dtemp_profile)  {
#
#     with(as.list(c(parms, x)), {
#
#         temp <- temp_profile(t)
#         dtemp <- dtemp_profile(t)
#
#         D_T <- D_R * 10^( (temp_ref - temp)/z )
#
#         correction <- abs( log(10)/z * t * dtemp + 1 )
#
#         dN <- - N * log(10)/D_T * correction
#
#         res <- c(dN)
#         return(list(res))
#
#     })
# }

