
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA 02111-1307 USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  maxddStats                Expectation of Drawdowns for BM with drift
#  .maxddStats               Utility function called by "maxddStats"
# FUNCTION:                 DISTRIBUTION AND RANDOM VARIATES:
#  dmaxdd                    Density function of mean Max-Drawdowns
#  pmaxdd                    Probability function of mean Max-Drawdowns
#  rmaxdd                    Random Variates of mean Max-Drawdowns
################################################################################


maxddStats <- 
function(mean = 0, sd = 1, horizon = 1000)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates Expectation Value E[D] of maximum Drawdowns of
    #   Brownian Motion for a given drift "mu", variance "sigma",  
    #   and runtime "horizon" of the Brownian Motion process.
    
    # Arguments:
    #   mu - Drift of Brownian Motion, a numeric vector.
    #   sigma - Standard Deviation of Brownian Motion, a numeric value.
    #   horizon - Runtime of the process, the time horizon of the 
    #     investor, a numeric value.
    
    # Details:
    #   Interpolates the data given in the table of Appendix B
    #   in Magdon-Ismail et al. [2003], "On the Maximum Drawdown 
    #   of Brownian Motion". The interpolation is done with R's
    #   function spline().
    
    # Notes.
    #   This computes for a vector of horizon values.
    
    # FUNCTION:
    
    # Statistics:
    ans = NULL
    for (i in 1:length(horizon)) {
        ans = c(ans, .maxddStats(mu = mean, sigma = sd, horizon = horizon[i]))
    }
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.maxddStats <- 
function(mu = 0, sigma = 1, horizon = 1000)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Utility function called by "maxddStats"
    
    # Arguments:
    #   see function "maxddStats
    
    # FUNCTION:
    
    # Internal Function - POSITIVE CASE: mu > 0
    # Left Table from Appendix B:
    QP = function(x) {
        gamma = sqrt(pi/8)
        vqn = t(matrix(c(
            0.0005,0.019690, 0.0010,0.027694, 0.0015,0.033789, 0.0020,0.038896,
            0.0025,0.043372, 0.0050,0.060721, 0.0075,0.073808, 0.0100,0.084693,
            0.0125,0.094171, 0.0150,0.102651, 0.0175,0.110375, 0.0200,0.117503,
            0.0225,0.124142, 0.0250,0.130374, 0.0275,0.136259, 0.0300,0.141842,
            0.0325,0.147162, 0.0350,0.152249, 0.0375,0.157127, 0.0400,0.161817,
            0.0425,0.166337, 0.0450,0.170702, 0.0500,0.179015, 0.0600,0.194248,
            0.0700,0.207999, 0.0800,0.220581, 0.0900,0.232212, 0.1000,0.243050,
            0.2000,0.325071, 0.3000,0.382016, 0.4000,0.426452, 0.5000,0.463159,
            1.5000,0.668992, 2.5000,0.775976, 3.5000,0.849298, 4.5000,0.905305,
            10.000,1.088998, 20.000,1.253794, 30.000,1.351794, 40.000,1.421860,
            50.000,1.476457, 150.00,1.747485, 250.00,1.874323, 350.00,1.958037,
            450.00,2.020630, 1000.0,2.219765, 2000.0,2.392826, 3000.0,2.494109,
            4000.0,2.565985, 5000.0,2.621743), 2))
        # Interpolation:
        if (x < 0.0005) { y = gamma*sqrt(2*x) }
        if (x >= 0.0005 & x <= 5000) {
            y = spline(log(vqn[,1]), vqn[,2], n = 1, xmin = log(x), 
                xmax = log(x))$y 
        }
        if (x > 5000) { y = 0.25*log(x) + 0.49088}
        # Return Value:
        y 
    }
    
    # Internal Function - NEGATIVE CASE: mu < 0
    # Right Table from Appendix B:
    QN = function(x) {
        gamma = sqrt(pi/8)
        vqn = t(matrix(c(
            0.0005,0.019965, 0.0010,0.028394, 0.0015,0.034874, 0.0020,0.040369,
            0.0025,0.045256, 0.0050,0.064633, 0.0075,0.079746, 0.0100,0.092708,
            0.0125,0.104259, 0.0150,0.114814, 0.0175,0.124608, 0.0200,0.133772,
            0.0225,0.142429, 0.0250,0.150739, 0.0275,0.158565, 0.0300,0.166229,
            0.0325,0.173756, 0.0350,0.180793, 0.0375,0.187739, 0.0400,0.194489,
            0.0425,0.201094, 0.0450,0.207572, 0.0475,0.213877, 0.0500,0.220056,
            0.0550,0.231797, 0.0600,0.243374, 0.0650,0.254585, 0.0700,0.265472,
            0.0750,0.276070, 0.0800,0.286406, 0.0850,0.296507, 0.0900,0.306393,
            0.0950,0.316066, 0.1000,0.325586, 0.1500,0.413136, 0.2000,0.491599,
            0.2500,0.564333, 0.3000,0.633007, 0.3500,0.698849, 0.4000,0.762455,
            0.5000,0.884593, 1.0000,1.445520, 1.5000,1.970740, 2.0000,2.483960,
            2.5000,2.990940, 3.0000,3.492520, 3.5000,3.995190, 4.0000,4.492380,
            4.5000,4.990430, 5.0000,5.498820), 2))
        # Interpolation:
        if (x < 0.0005) { y = gamma*sqrt(2*x) }
        if (x >= 0.0005 & x <= 5) {
            y = spline(vqn[,1], vqn[,2], n = 1, xmin = x, xmax = x)$y }
        if (x > 5) { y = x + 1/2 }  
        # Return Value:
        y 
    }
    
    # Result:
    ED = NULL
    for (i in 1:length(mu)) {
        if (mu[i] == 0) {
            gamma = sqrt(pi/8)  
            ED[i] = 2 * gamma * sigma * sqrt(horizon) 
        } else {
            x = mu[i]^2 * horizon / ( 2 * sigma^2 )
            if (mu[i] > 0) { ED[i] = ( 2* sigma^2 / mu[i] ) * QP(x) }
            if (mu[i] < 0) { ED[i] = - ( 2* sigma^2 / mu[i] ) * QN(x) } 
        } 
    }
    
    # Return Value:
    ED
}
    

################################################################################

    
dmaxdd <- 
function(x, sd = 1, horizon = 100, N = 1000) 
{
    # Description:
    #   Calculates for a trendless Brownian process (mean=0) and
    #   standard deviation "sd" the density from the probability 
    #   that the maximum drawdown "D" is larger or equal to "h"  
    #   in the interval [0,T], where "T" denotes the time horizon  
    #   of the investor.
    
    # Arguments:
    #   x - Vector of Drawdowns
    #   sd - Standard Deviation
    #   horizon - Time horizon of the investor
    #   N - number of summands
    
    # Notes:
    #   The drift dependent case is not yet implemented!
    
    # FUNCTION:
        
    # Use "h", "sigma" to be conform with Magdon-Ismael et al. [2003]
    h = x
    sigma = sd
    
    # Settings:
    n = 1:N
    pn = 2 * sin((n-0.5)*pi) * sigma^2 * (n-0.5) * pi * horizon
    en = sigma^2 * (n-0.5)^2 * pi^2 * horizon / 2
    
    # Loop over all Drawdowns:
    result = rep(NA, times = length(h))         
    for (i in 1:length(h)) {
        if (h[i] == 0) {
            result[i] = 0 
        } else {
            g = pn *  exp(-en/(h[i]^2)) / h[i]^3
            result[i] = sum(g) 
        } 
    }
    
    # Return Value:
    result
}   


# ------------------------------------------------------------------------------


pmaxdd <- 
function(q, sd = 1, horizon = 100, N = 1000) 
{
    # Description:
    #   Calculates for a trendless Brownian process (mean=0) 
    #   with standard deviation "sd" the probability that the 
    #   maximum drawdown "D" is larger or equal to "h" in the 
    #   interval [0,T], where "T" denotes the time horizon of 
    #   the investor.
    
    # Arguments:
    #   q - Vector of Drawdowns
    #   sd - Standard Deviation
    #   horizon - Time horizon of the investor
    #   N - number of summands
    
    # Value:
    #   GD(h) - eqn(14) In Magdon-Ismail et al.
    
    # Notes:
    #   The drift dependent case is not yet implemented!
    
    # FUNCTION:
        
    # Use "h", "sigma" to be conform with Magdon-Ismael et al. [2003]
    h = q
    sigma = sd
    
    # Settings:
    n = 1:N
    pn = 2 * ( sin((n-0.5)*pi) / ((n-0.5)*pi) ) 
    en = sigma^2 * (n-0.5)^2 * pi^2 * horizon / 2
    
    # Loop over all Drawdowns:
    result = rep(NA, times = length(h))         
    for (i in 1:length(h)) {
        g = pn * ( 1 - exp(-en/(h[i]^2)) )
        result[i] = sum(g) 
    }
    
    # Return Value:
    result
}   


# ------------------------------------------------------------------------------


rmaxdd <- 
function(n, mean = 0, sd = 1, horizon = 100) 
{
    # Description:
    #   Generates for a Brownian process with mean "mean" and
    #   standard deviation "sd" random variates of maximum
    #   Drawdowns.
    
    # Arguments:
    #   n - Number of Drawdown rvs
    #   mean - Drift
    #   sd - Standard Deviation
    #   horizon - Time horizon of the investor
    
    # FUNCTION:
        
    # Simulation of "n" max Drawdowns "h":
    result = NULL
    for (i in 1:n) {
        D = cumsum(rnorm(horizon, mean = mean, sd = sd))
        result[i] = max(cummax(D)-D) 
    }
    
    # Return Value:
    result
}   


################################################################################

