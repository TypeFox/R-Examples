
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
# FUNCTION:           DESCRIPTION:
#  HNGOption           Computes Option Price from the HN-GARCH Formula
#  HNGGreeks           Calculates one of the Greeks of the HN-GARCH Formula
#  HNGCharacteristics  Computes Option Price and all Greeks of HN-GARCH Model
################################################################################


test.HNGOption = 
function()
{
    # HNGOption - Computes Option Price from the HN-GARCH Formula
    
    # Define the Model Parameters for a Heston-Nandi Option:
    model = list(lambda = -0.5, omega = 2.3e-6, alpha = 2.9e-6, 
        beta = 0.85, gamma = 184.25) 
    S = X = 100
    Time.inDays = 252
    r.daily = 0.05/Time.inDays
    sigma.daily = sqrt((model$omega + model$alpha) /
        (1 - model$beta - model$alpha * model$gamma^2))
    data.frame(S, X, r.daily, sigma.daily)
    
    # HNGOption:
    
    # Compute HNG Call-Put and compare with GBS Call-Put:
    HNG = GBS = Diff = NULL
    for (TypeFlag in c("c", "p")) {
        HNG = c(HNG, HNGOption(TypeFlag, model = model, S = S, X = X, 
            Time.inDays = Time.inDays, r.daily = r.daily)$price )
        GBS = c(GBS, GBSOption(TypeFlag, S = S, X = X, Time = Time.inDays, 
            r = r.daily, b = r.daily, sigma = sigma.daily)@price) 
    }
    Options = cbind(HNG, GBS, Diff = round(100*(HNG-GBS)/GBS, digits = 2))
    row.names(Options) <- c("Call", "Put")
    data.frame(Options)
    
    # TODO: HNG not yet a S4 Class Member !!!

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.HNGGreeks = 
function()
{       
    # HNGGreeks - Calculates one of the Greeks of the HN-GARCH Formula
    
    # Define the Model Parameters for a Heston-Nandi Option:
    model = list(lambda = -0.5, omega = 2.3e-6, alpha = 2.9e-6, 
        beta = 0.85, gamma = 184.25) 
    S = X = 100
    Time.inDays = 252
    r.daily = 0.05/Time.inDays
    sigma.daily = sqrt((model$omega + model$alpha) /
        (1 - model$beta - model$alpha * model$gamma^2))
    data.frame(S, X, r.daily, sigma.daily)
    
    # Compute HNG Greeks and compare with GBS Greeks:
    Selection = c("Delta", "Gamma")
    HNG = GBS = NULL
    for (i in 1:2){
        HNG = c(HNG, HNGGreeks(Selection[i], TypeFlag = "c", model = model, 
            S = 100, X = 100, Time = Time.inDays, r = r.daily)) 
        GBS = c(GBS, GBSGreeks(Selection[i], TypeFlag = "c", S = 100, X = 100, 
            Time = Time.inDays, r = r.daily, b = r.daily, sigma = sigma.daily))
    }
    Greeks = cbind(HNG, GBS, Diff = round(100*(HNG-GBS)/GBS, digits = 2))
    row.names(Greeks) <- Selection
    data.frame(Greeks)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.HNGCharacteristics = 
function()
{ 
    # HNGCharacteristics  
    #   Computes Option Price and all Greeks of HN-GARCH Model
    
    NA
    
    # Return Value:
    return()    
}


################################################################################

