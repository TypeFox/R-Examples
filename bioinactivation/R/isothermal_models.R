
#' Isothermal Bigelow's Model
#'
#' Returns the predicted logarithmic reduction in microbial count according
#' to Bigelow's model for the time, temperature and model parameters given.
#'
#' @param time numeric indicating the time at which the prediction is taken.
#' @param temp numeric indicating the temperature of the treatment.
#' @param z numeric defining the z-value.
#' @param D_R numeric defining the D-value at the reference temperature.
#' @param temp_ref numeric defining the reference temperature.
#'
#' @return A numeric with the predicted logarithmic reduction
#'         (\eqn{log10(N/N0)}).
#'
Bigelow_iso <- function(time, temp, z, D_R, temp_ref){

    D_T <- D_R * 10^( (temp_ref - temp)/z )
    log_diff <- -time/D_T
    return(log_diff)
}

#' Isothermal Weibull-Mafart Model
#'
#' Returns the predicted logarithmic reduction in microbial count according
#' to Weibull-Mafarts's model for the time, temperature and model parameters
#' given.
#'
#' @param time numeric indicating the time at which the prediction is taken.
#' @param temp numeric indicating the temperature of the treatment.
#' @param delta_ref numeric defining the delta-value at the reference temperature.
#' @param z numeric defining the z-value.
#' @param p numeric defining shape factor of the Weibull distribution.
#' @param temp_ref numeric indicating the reference temperature.
#'
#' @return A numeric with the predicted logarithmic reduction
#'         (\eqn{log10(N/N0)}).
#'
WeibullMafart_iso <- function(time, temp, delta_ref, z, p, temp_ref){

    delta <- delta_ref * 10^(- (temp - temp_ref)/z )
    log_diff <- -(time/delta)^p
    return(log_diff)
}

#' Isothermal Weibull-Peleg Model
#'
#' Returns the predicted logarithmic reduction in microbial count according
#' to Weibull-Peleg's model for the time, temperature and model parameters
#' given.
#'
#' @param time numeric indicating the time at which the prediction is taken.
#' @param temp numeric indicating the temperature of the treatment.
#' @param n numeric defining shape factor of the Weibull distribution.
#' @param k_b numeric indicating the slope of the b~temp line.
#' @param temp_crit numeric with the value of the critical temperature.
#'
#' @return A numeric with the predicted logarithmic reduction
#'         (\eqn{log10(N/N0)}).
#'
WeibullPeleg_iso <- function(time, temp, n, k_b, temp_crit){

    b <- log( 1 + exp(k_b*(temp - temp_crit)))
    log_diff <- -b * time^n
    return(log_diff)

}

### ---------------------------------------------------------- ###
### SECONDARY  MODEL FOR N_RES IS NEEDED  TO DEVELOP THE MODEL ###

# #' Isothermal Geeraerd Model
# #'
# #' Returns the predicted logarithmic reduction in microbial cont according
# #' to Geeraerd's model.
# #' The isothermal prediction is calculated by analytical integration of the
# #' ode for constant temperature
# #'
# Geeraerd_iso <- function(time, temp, C_c0, N_min, z, D_R, temp_ref) {
#
#     D_T <- D_R * 10^( (temp_ref - temp)/z )
#     k_max <- 1/D_T
#     N0 <- 1e6
#     C1 <- (C_c0 + 1)*N0 - N_min
#     N <- C1/(C_c0+exp(k_max*time)) + N_min*exp(k_max*time)/(C_c0+exp(k_max*time))
#     log_diff <- log10(N/N0)
#     return(log_diff)
#
#
# }

