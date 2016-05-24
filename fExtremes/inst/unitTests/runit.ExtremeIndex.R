
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
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPTION:
#  'fTHETA'                  Class representation for extremal index
#  show.fTHETA               S4: Print Method for extremal index
#  thetaSim                  Simulates a time series with known theta
# FUNCTION:                 DESCRIPTION:
#  blockTheta                Computes theta from Block Method
#  clusterTheta              Computes theta from Reciprocal Cluster Method
#  runTheta                  Computes theta from Run Method
#  ferrosegersTheta          Computes Theta according to Ferro and Seegers
# FUNCTION:                 DESCRIPTION:
#  exindexesPlot             Computes and Plot Theta(1,2,3)
#  exindexPlot               Computes Theta(1,2) and Plot Theta(1)
################################################################################


test.fTHETA = 
function()
{
    # Slot Names:
    slotNames("fTHETA")
    # [1] "call" "data" "theta" "title" "description"

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.thetaSim = 
function()
{
    # Simulation:
    # thetaSim(model = c("max", "pair"), n = 100, theta = 0.5)

    # Max Frechet Series:
    x = thetaSim("max")
    class(x)
    print(x)
    
    # Paired Exponential Series:
    x = thetaSim("pair")
    class(x)
    print(x)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.thetaFit = 
function()
{
    # Parameter Estimation:
    x.ts = thetaSim("max", n=22000)
    class(x.ts)
    
    # Parameter Estimation:
    # blockTheta(x, block = 22, quantiles = seq(0.95, 0.995, length = 10), 
    #   title = NULL, description = NULL) 
    # clusterTheta(x, block = 22, quantiles = seq(0.95, 0.995, length = 10), 
    #   title = NULL, description = NULL) 
    # runTheta(x, block = 22, quantiles = seq(0.95, 0.995, length = 10), 
    #   title = NULL, description = NULL)  
    # ferrosegersTheta(x, quantiles = seq(0.95, 0.995, length = 10), 
    #   title = NULL, description = NULL) 
       
    # time series ts as input:
    blockTheta(x.ts)
    clusterTheta(x.ts)
    runTheta(x.ts)
    ferrosegersTheta(x.ts)
    
    # Numeric Vector as input:
    x.vec = as.vector(x.ts)
    blockTheta(x.vec)
    clusterTheta(x.vec)
    runTheta(x.vec)
    ferrosegersTheta(x.vec)
    
    # timeSeries object as input:
    x.tS = as.timeSeries(x.ts)
    blockTheta(x.tS)
    clusterTheta(x.tS)
    runTheta(x.tS)
    ferrosegersTheta(x.tS)      
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.exindexesPlot = 
function()
{
    # Graphics Frame:
    par(mfrow = c(2, 2), cex = 0.7)
    par(ask = FALSE)
    
    # Parameter Estimation:
    x = thetaSim("max", n = 22000)
    exindexesPlot(x)
    
    # Parameter Estimation:
    y = thetaSim("pair", n = 22000)
    exindexesPlot(y)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.exindexPlot = 
function()
{
    # Graphics Frame:
    par(mfrow = c(2, 2), cex = 0.7)
    par(ask = FALSE)
    
    # Parameter Estimation:
    x = thetaSim("max", n=22000)
    exindexPlot(x, block = 22)
    
    # Parameter Estimation:
    y = thetaSim("pair", n=22000)
    exindexPlot(y, block = 22)
    
    # Return Value:
    return()    
}


################################################################################

