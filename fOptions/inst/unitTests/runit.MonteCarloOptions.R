
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                  DESCRIPTION:
#  MonteCarloOption           Valuate Options by Monte Carlo Simulation
################################################################################


test.MonteCarloOption <- 
    function()
{
    # How to perform a Monte Carlo Simulation?
    
    # RVs:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    
    # First Step:
    # Write a function to generate the option's innovations. 
    # Use scrambled normal Sobol numbers:
    sobolInnovations = function(mcSteps, pathLength, init, ...) {
        # Create Normal Sobol Innovations:
        innovations = rnorm.sobol(mcSteps, pathLength, init, ...)
        # Return Value:
        innovations 
    }
    
    # Second Step: 
    # Write a function to generate the option's price paths.  
    # Use a Wiener path:
    wienerPath = function(eps) { 
        # Note, the option parameters must be globally defined!
        # Generate the Paths:
        path = (b-sigma*sigma/2)*delta.t + sigma*sqrt(delta.t)*eps
        # Return Value:
        path 
    }
      
    # Third Step: 
    # Write a function for the option's payoff
    
    # Example 1: use the payoff for a plain Vanilla Call or Put:
    plainVanillaPayoff = function(path) { 
        # Note, the option parameters must be globally defined!
        # Compute the Call/Put Payoff Value:
        ST = S*exp(sum(path))
        if (TypeFlag == "c") payoff = exp(-r*Time)*max(ST-X, 0)
        if (TypeFlag == "p") payoff = exp(-r*Time)*max(0, X-ST)
        # Return Value:
        payoff 
    }
    
    # Example 2: use the payoff for an arithmetic Asian Call or Put:
    arithmeticAsianPayoff = function(path) { 
        # Note, the option parameters must be globally defined!
        # Compute the Call/Put Payoff Value:
        SM = mean(S*exp(cumsum(path)))
        if (TypeFlag == "c") payoff = exp(-r*Time)*max(SM-X, 0)
        if (TypeFlag == "p") payoff = exp(-r*Time)*max(0, X-SM)
        # Return Value:
        payoff 
    }
    
    # Final Step: 
    # Set Global Parameters for the plain Vanilla / arithmetic Asian Options:
    TypeFlag <- "c"; S <- 100; X <- 100
    Time <- 1/12; sigma <- 0.4; r <- 0.10; b <- 0.1
    
    # Do the Asian Simulation with scrambled random numbers:
    mc = MonteCarloOption(delta.t = 1/360, pathLength = 30, mcSteps = 5000, 
        mcLoops = 50, init = TRUE, innovations.gen = sobolInnovations, 
        path.gen = wienerPath, payoff.calc = arithmeticAsianPayoff, 
        antithetic = TRUE, standardization = FALSE, trace = TRUE, 
        scrambling = 2, seed = 4711)
    
    # Plot the MC Iteration Path:
    par(mfrow = c(1, 1))
    mcPrice = cumsum(mc)/(1:length(mc))
    plot(mcPrice, type = "l", main = "Arithmetic Asian Option", 
        xlab = "Monte Carlo Loops", ylab = "Option Price")
    
    # Compare with Turnbull-Wakeman Approximation:
    #   ... requires(fExoticOptions)
    # TW = TurnbullWakemanAsianApproxOption(TypeFlag = "c", S = 100, SA = 100, 
    #    X = 100, Time = 1/12, time = 1/12, tau = 0 , r = 0.1, b = 0.1, 
    #    sigma = 0.4)$price
    # print(TW)
    TW = 2.859122 
    abline(h = TW, col = 2)

    # Return Value:
    return()    
}


################################################################################

