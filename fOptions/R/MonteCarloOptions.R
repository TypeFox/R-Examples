
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


################################################################################
# FUNCTION:                  DESCRIPTION:
#  MonteCarloOption           Valuate Options by Monte Carlo Simulation
################################################################################


MonteCarloOption = function(delta.t, pathLength, mcSteps, mcLoops, 
init = TRUE, innovations.gen, path.gen, payoff.calc, antithetic = TRUE, 
standardization = FALSE, trace = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Valuates Options by Monte Carlo Simulation

    # Arguments:
    #   delta.t    - The length of the time interval, by default one day.
    #   pathLength - Number of Time Intervals which add up to the path.
    #   mcSteps    - The number of Monte Carlo Steps performed in one loop.
    #   mcLoops    - The number of Monte Carlo Loops
    #   init       - Should the random number generator be initialized ?
    #                This variable must appear in the argument list of the
    #                random number generator, even it will not ne used
    #   innovations.gen 
    #              - the generator function for the innovations
    #   path.gen   - the generator for the MC paths
    #   payoff.calc
    #              - the payoff claculator function
    #   antithetic - if TRUE, antithetic paths are used in the MC simulation
    #   standardization
    #              - if TRUE, the random numbers will be standardized so that
    #                their mean is zero and their variance is zero
    #   trace      - a logical, should the iteration path be traced ?
    #   ...        - additional parameters passed to innovations.gen.
    
    # Notes:
    #   Global Variables:
    #     The options parameter must be globally available. 
    #     For a Black-Scholes or a simple Asian Option these are:
    #     TypeFlag, S, X, Time, r, b, sigma
    #   Required Functions:
    #     The user must specify the following functions:
    #     innovations.gen()
    #     path.gen()
    #     payoff.calc()
    
    # FUNCTION
        
    # Monte Carlo Simulation:
    delta.t <<- delta.t
    if (trace) cat("\nMonte Carlo Simulation Path:\n\n")
        iteration = rep(0, length = mcLoops)
        # MC Iteration Loop:
        cat("\nLoop:\t", "No\t")
        for ( i in 1:mcLoops ) {
            if ( i > 1) init = FALSE
            # Generate Innovations:
              eps =  innovations.gen(mcSteps, pathLength, init = init, ...)
            # Use Antithetic Variates if requested:
              if (antithetic) 
                  eps = rbind(eps, -eps)
            # Standardize Variates if requested:
              if (standardization) eps = 
                 (eps-mean(eps))/sqrt(var(as.vector(eps)))
            # Calculate for each path the option price:
              path = t(path.gen(eps))
              payoff = NULL
              for (j in 1:dim(path)[1]) 
                  payoff = c(payoff, payoff.calc(path[, j])) 
              iteration[i] = mean(payoff)
            # Trace the Simualtion if desired:
              if (trace) 
                 cat("\nLoop:\t", i, "\t:", iteration[i], sum(iteration)/i ) 
        }
        if (trace) cat("\n")

    # Return Value:
    iteration
}


# ******************************************************************************

