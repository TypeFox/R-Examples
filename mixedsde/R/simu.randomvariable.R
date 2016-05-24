#' Simulation Of Random Variables
#' 
#' @description Simulation of (discrete) random variables from a vector of probability (the nonparametrically estimated values
#'  of the density renormalised to sum at 1) and a vectors of real values (the grid of estimation)
#' @param x n real numbers
#' @param p vector of probability, length n

#' @return
#' y a simulated value from the discrete distribution



discr <- function(x, p) {
    cp <- cumsum(p)
    U <- runif(1, 0, max(cp))
    x[which(cp >= U)[1]]
} 
