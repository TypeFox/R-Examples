
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
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################


test.plot.methods1 <- 
    function()
{
    # Load data:
    data(dem2gbp)
    dem2gbp = as.vector(dem2gbp[, 1])
    
    # Fit to normal Conditional Distribution:
    fit = garchFit( ~ garch(1, 1), data = dem2gbp, trace = FALSE)
    print(fit)
    
    # garchFit 1:
    # 1:
    # 2:
    
    # Graph Frame:
    par(mfrow = c(2, 1))
    
    # Plot 1:
    plot(fit, which = 1)                      
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, col = "darkgrey")
    
    # Plot 2:
    plot(fit, which = 2)
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, col = "darkgrey")
    
   
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.plot.methods2 <- 
    function()
{
    # Load data:
    data(dem2gbp)
    
    # Fit to normal Conditional Distribution:
    fit = garchFit( ~ garch(1, 1), data = dem2gbp, trace = FALSE)
    print(fit)
    
    # garchFit 2:
    
    # Graph Frame:
    par(mfrow = c(1, 1))
    
    # Plot 3:
    plot(fit, which = 3)                       
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")

    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.plot.methods3 <- 
    function()
{
    # Load data:
    data(dem2gbp)
    
    # Fit to normal Conditional Distribution:
    fit = garchFit( ~ garch(1, 1), data = dem2gbp, trace = FALSE)
    print(fit)
    
    # garchFit3:
    # 3:
    # 4:
    # 5:
    
    # Graph Frame:
    par(mfrow = c(2, 1))
    
    # Plot 4:
    plot(fit, which = 4)                       
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")
    
    # Plot 5:
    plot(fit, which = 5)
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")

    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.plot.methods4 <- 
    function()
{
    # Load data:
    data(dem2gbp)
    
    # Fit to normal Conditional Distribution:
    fit = garchFit( ~ garch(1, 1), data = dem2gbp, trace = FALSE)
    print(fit)
    
    # garchFit4:
    # 6: Cross Correlation
    # 7: Residuals
    # 8: Conditional SDs
    # 9: Standardized Residuals
    
    # Graph Frame:
    par(mfrow = c(2, 2))
    
    # Plot 6:
    plot(fit, which = 6)                        
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")
    
    # Plot 7:
    plot(fit, which = 7)
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")
    
    # Plot 8:
    plot(fit, which = 8)                       
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")
    
    # Plot 9:
    plot(fit, which = 9)
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")
   
    # Return Value:
    return()    
} 


# ------------------------------------------------------------------------------


test.plot.methods5 <- 
    function()
{
    # Load data:
    data(dem2gbp)
    
    # Fit to normal Conditional Distribution:
    fit = garchFit( ~ garch(1, 1), data = dem2gbp, trace = FALSE)
    print(fit)
    
    # garchFit5:
    # 10: ACF of Standardized Residuals
    # 11: ACF of Squared Standardized Residuals
    # 12: Cross Correlation between r^2 and r
    # 13: QQ-Plot of Standardized Residuals
    
    # Graph Frame:
    par(mfrow = c(2, 2))
    
    # Plot 10:
    plot(fit, which = 10)                      
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")
    
    # Plot 11:
    plot(fit, which = 11)
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")
    
    # Plot 12:
    plot(fit, which = 12)                       
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")
    
    # Plot 13:
    plot(fit, which = 13)
    mtext("norm-GARCH(1,1) Modeling", line = 0.5, cex = 0.8)
    mtext("DEM2GBP Data Vector", side = 4, adj = 0, cex = 0.7, 
        col = "darkgrey")

    # Return Value:
    return()    
} 


################################################################################

