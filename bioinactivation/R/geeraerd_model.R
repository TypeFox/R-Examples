
#'
#' First Derivate of Geeraerd's Model
#'
#' Calculates the first derivative of the Geeraerd's model at a given time for
#' the model parameters provided and the environmental conditions given.
#'
#' This function is compatible with the function
#' \code{\link{predict_inactivation}}.
#'
#' @section Model Equation:
#'      \deqn{\frac{dN}{dt} = -N \cdot k_{max} \cdot \alpha \cdot \gamma}{
#'            dN/dt = - N * k_max * \alpha * \gamma}
#'
#'      \deqn{\frac{dC_c}{dt} = - C_c \cdot k_{max}}{
#'            dC_c/dt = - C_c * k_max}
#'
#'      \deqn{\alpha = \frac {1}{1+C_c}}{
#'            \alpha = 1/(1+C_c)}
#'
#'      \deqn{\gamma = 1 - \frac{N_{res}}{N}}{
#'            \gamma = 1 - N_res/N}
#'
#'      \deqn{k_{max} = \frac{1}{D_T}}{k_max = 1/D_T}
#'
#'      \deqn{log_{10}D_T = log_{10}D_R + \frac{T_R-T}{z}}{
#'             log(D_T) = log(D_R) + (T_R-T)/z}
#'
#' @param t numeric vector indicating the time of the experiment.
#' @param x list with the values of the variables at time \eqn{t}.
#' @param parms parameters for the secondary model. No explicit check of their
#'        validity is performed (see section \bold{Model Parameters}).
#' @param temp_profile a function that provides the temperature at a given time.
#'
#' @return A list with the value of the first derivatives of N and C_c at time
#'         \code{t}.
#'
#' @section Model Parameters:
#'      \itemize{
#'          \item temp_ref: Reference temperature for the calculation,
#'          \item D_R: D-value at the reference temperature,
#'          \item z: z value,
#'          \item N_min: Minimum value of N (defines the tail).
#'          }
#'
#' @section Notes:
#'      To define the Geeraerd model without tail, assign \code{N_min = 0}.
#'      For the model without shoulder, assign \code{C_0 = 0}
#'
#' @seealso \code{\link{predict_inactivation}}
#'
dGeeraerd_model<- function(t, x, parms, temp_profile)  {

    temp <- temp_profile(t)

    with(as.list(c(x, parms)),{

        D_T <- D_R*10^((temp_ref - temp)/z)
        k_max <- log(10)/D_T

        dC_c <- - k_max * C_c
        dN <- - N * k_max * (1/(1+C_c)) * (1 - N_min/N)

        res <- c(dN, dC_c)
        return(list(res))
    })

}

