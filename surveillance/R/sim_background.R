###################################################
### chunk number 1: 
###################################################

# 'sim.seasonalNoise' generates a cyclic model of a poisson distribution
# as background data for a simulated timevector.
#
# Parameters:
# A             - amplitude (range of sinus), default = 1
# alpha         - parameter to move along the y-axis
#                       (negative values not allowed)
#                       d.h alpha > = A, default = 1,
# beta          - regression coefficient, default = 0
# season        - factor to create seasonal moves
#                       (moves the curve along the x-axis), default = 0
# length                - number of weeks to model
# frequency     - factor to determine the oscillation-frequency, default = 1
# state        - if a state chain is given, it is weighted by the parameter K
#                       and influences mu
# K             - weight for outbreaks

sim.seasonalNoise <- function(A = 1, alpha = 1, beta = 0, phi = 0,
                                length, frequency = 1, state = NULL, K = 0){

        t <- 1:length

        # constant factor to transform weeks to the appropriate pi-value.
        omega <- 2 * pi/ 52

        # season moves the sin along the x-axis.
        if(is.null(state)){ # no state chain
                mu <- exp(A * sin( frequency * omega * (t + phi)) + alpha  + beta * t)
        }
        else{ # encounter the state chain
                mu <- exp(A * sin( frequency * omega * (t + phi)) +
                                        alpha  + beta * t +  K * state)
        }

        # create the noise as random numbers of the Poisson distribution
        # with parameter mu
        seasonalBackground <- rpois(length, mu) # get random numbers

        result <- list(seasonalBackground = seasonalBackground, t = t,
                                mu = mu, A = A, alpha = alpha,
                                beta = beta, phi = phi,
                                length = length,
                                frequency = frequency, K = K)

        class(result) = "seasonNoise"
        return(result)
}






