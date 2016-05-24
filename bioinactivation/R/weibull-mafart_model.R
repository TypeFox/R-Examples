
#'
#' First Derivate of the Weibull-Mafart Model
#'
#' Calculates the first derivative of Weibull-Mafart model at a given time for
#' the model parameters provided and the environmental conditions given.
#'
#' The model is developed from the isothermal Weibull-Mafart model without
#' taking into
#' account in the derivation the time dependence of \eqn{\delta_T} for
#' non-isothermal temperature profiles.
#'
#' This function is compatible with the function
#' \code{\link{predict_inactivation}}.
#'
#' @section Model Equation:
#'
#'      \deqn{\frac{dN}{dt} = -N \cdot p \cdot (1/\delta)^p \cdot t^{p-1} }{
#'            dN/dt = -N * p * (1/delta)^p * t^(p-1)}
#'
#'      \deqn{\delta(T) = \delta_{ref} \cdot 10^{- (T-T_ref)/z} }{
#'            delta(T) = delta_ref * 10^(- (T-T_ref)/z )}
#'
#' @param t numeric vector indicating the time of the experiment.
#' @param x list with the value of N at t.
#' @param parms parameters for the secondary model. No explicit check of their validity
#'             is performed (see section \bold{Model Parameters}).
#' @param temp_profile a function that provides the temperature at a given time.
#'
#' @return The value of the first derivative of N at time \code{t} as a list.
#'
#' @section Model Parameters:
#'      \itemize{
#'          \item temp_ref: Reference temperature for the calculation.
#'          \item delta_ref: Value of the scale factor at the reference temperature.
#'          \item z: z-value.
#'          \item p: shape factor of the Weibull distribution.
#'          }
#'
#' @section Note:
#'      For t=0, dN = 0 unless n=1. Hence, a small shift needs to be introduced
#'      to t.
#'
#' @seealso \code{\link{predict_inactivation}}
#'
dMafart_model<- function(t, x, parms, temp_profile)  {

    temp <- temp_profile(t)

    with(as.list(c(x, parms)),{
        delta <- delta_ref * 10^( -(temp-temp_ref)/z)
        dN <- - N * p * (1/delta)^p * t^(p-1) * log(10)

        res <- c(dN)
        return(list(res))
    })

}

# #'
# #' First Derivate of the Mafart Model with Full Derivation
# #'
# #' Calculates the first derivative of Weibull-Mafart model at a given time for
# #' the model parameters provided and the environmental conditions given.
# #'
# #' The model is developed from the isothermal Weibull-Mafart model taking into
# #' account in the derivation the time dependence of \eqn{\delta_T} for
# #' non-isothermal temperature profiles.
# #' For the linearized version, see \code{\link{dMafart_model}}
# #'
# #' This function is compatible with the function
# #' \code{\link{predict_inactivation}}.
# #'
# #' @section Model Equation:
# #'
# #'      \deqn{ \frac{dN}{dt} = -N \cdot p \cdot (t/\delta)^{p-1} \cdot
# #'                \frac{\mathrm{ln}(10)}{\delta} (1+t \frac{\mathrm{ln}(10)}{z}
# #'                    \frac{dT}{dt}) }{
# #'            dN/dt =-N*p*(t/delta)^(p-1)*log(10)/delta*(1+log(10)*t/z*dtemp)}
# #'
# #'      \deqn{\delta(T) = \delta_{ref} \cdot 10^{- (T-T_ref)/z}}{
# #'            delta(T) = delta_ref * 10^(- (T-T_ref)/z )}
# #'
# #' @param t numeric vector indicating the time of the experiment.
# #' @param x list with the value of N at t.
# #' @param parms parameters for the secondary model. No explicit check of their validity
# #'             is performed (see section \bold{Model Parameters}).
# #' @param temp_profile a function that provides the temperature at a given time.
# #' @param dtemp_profile a function that provides the first derivative of the
# #'        temperature at a given time.
# #'
# #' @return The value of the first derivative of N at time \code{t} as a list.
# #'
# #' @section Model Parameters:
# #'      \itemize{
# #'          \item temp_ref: Reference temperature for the calculation.
# #'          \item delta_ref: Value of the scale factor at the reference temperature.
# #'          \item z: z-value.
# #'          \item p: shape factor of the Weibull distribution.
# #'          }
# #'
# #' @section Note:
# #'      For t=0, dN = 0 unless n=1. Hence, a small shift needs to be introduced
# #'      to t.
# #'
# #' @seealso \code{\link{predict_inactivation}}
# #'
# dMafart_model_full<- function(t, x, parms, temp_profile, dtemp_profile)  {
#
#     temp <- temp_profile(t)
#     dtemp <- dtemp_profile(t)
#
#     with(as.list(c(x, parms)),{
#         delta <- delta_ref * 10^( -(temp-temp_ref)/z)
#         dN <- -N * p * (t/delta)^(p-1) * log(10)/delta*abs( 1+log(10)*t/z*dtemp )
#         # dN <- -N * p * (t/delta)^(p-1) * log(10)/delta*( 1+log(10)*t/z*abs(dtemp) )
#
#         res <- c(dN)
#         return(list(res))
#     })
#
# }

