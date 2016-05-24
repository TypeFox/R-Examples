
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
# FUNCTION:             CHAOTIC TIME SERIES MAPS:
#  tentSim               Simulates series from Tent map
#  henonSim              Simulates series from Henon map 
#  ikedaSim              Simulates series from Ikeda map
#  logisticSim           Simulates series from Logistic map
#  lorentzSim            Simulates series from Lorentz map
#  roesslerSim           Simulates series from Roessler map
# FUNCTION:             PHASE SPACE REPRESENTATION:
#  mutualPlot            Creates mutual information plot
#  fnnPlot               Creates false nearest neigbours plot
# FUNCTION:             NON STATIONARITY PLOTS:
#  recurrencePlot        Creates recurrence plot
#  separationPlot        Creates space-time separation plot
# FUNCTION:             LYAPUNOV EXPONENTS:
#  lyapunovPlot          Maximum Lyapunov plot               
################################################################################


test.tentSim = 
function()
{  
    #  tentSim - Simulates series from Tent map 
    
    # Tent Map:  
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    par (mfrow = c(1, 1))
    ts = tentSim(n = 1000, n.skip = 100, parms = c(a = 2), start = runif(1), 
        doplot = TRUE) 
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.henonSim = 
function()
{  
    #  henonSim - Simulates series from Henon map 
     
    # Henon Map - 2D:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    par (mfrow = c(1, 1))
    ts = henonSim(n = 1000, n.skip = 100, parms = c(a = 1.4, b = 0.3), 
        start = runif(2), doplot = TRUE) 
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ikedaSim = 
function()
{   
    #  ikedaSim - Simulates series from Ikeda map 

    # Ikeda Map - 2D:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    par (mfrow = c(2, 2))
    ts = ikedaSim(n = 1000, n.skip = 100, parms = c(a = 0.4, b = 6, c = 0.9), 
        start = runif(2), doplot = TRUE) 
    head(ts)
  
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.logisticSim = 
function()
{  
    #  logisticSim - Simulates series from Logistic map
    #  lorentzSim - Simulates series from Lorentz map
    #  roesslerSim - Simulates series from Roessler map  
    
    # Logistic Map:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    par (mfrow = c(1, 1))
    logisticSim(n = 1000, n.skip = 100, parms = c(r = 4), start = runif(1), 
        doplot = TRUE) 
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.lorentzSim = 
function()
{  
    #  lorentzSim - Simulates series from Lorentz map
    
    # Lorentz Map:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    par (mfrow = c(3, 2))
    ts = lorentzSim(times = seq(0, 20, by = 0.01), parms = c(sigma = 16, 
        r = 45.92, b = 4), start = c(-14, -13, 47), doplot = TRUE) 
    head(ts)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.roesslerSim = 
function()
{  
    #  roesslerSim - Simulates series from Roessler map  
    
    # Roessler Map:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    par (mfrow = c(3, 2))
    ts = roesslerSim(times = seq(0, 80, by = 0.05), parms = c(a = 0.2, 
        b = 0.2, c = 8), start = c(-1.894, -9.92, 0.025), doplot = TRUE) 
    head(ts)
    
    # Return Value:
    return()    
}


################################################################################


test.henonSlider = 
function()
{
    # Henon Slider:
    henonSlider = function()
    {  
        refresh.code = function(...)
        {
            # Sliders:
            N = .sliderMenu(no = 1)
            a = .sliderMenu(no = 2)
            b = .sliderMenu(no = 3)
            
            # Plot Henon Map:      
            ts = henonSim(n = N, n.skip = 100, parms = c(a = a, b = b), 
                start = c(pi/4, exp(1)/4), doplot = TRUE) 
            
            # Frame:
            par(mfrow = c(1, 1), cex = 0.7)
        }
      
        # Open Slider Menu:
        .sliderMenu(refresh.code,
           names =       c( "N",    "a",     "b"),
           minima =      c( 100,   1.00,    0.00),
           maxima =      c(5000,   2.00,    1.00),
           resolutions = c( 100,   0.01,    0.01),
           starts =      c(2000,   1.40,    0.30))
    }
    
    # Try:
    # henonSlider()
    
    # Return Value:
    return() 
}


################################################################################
    
