
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
# FUNCTION:                 HEAVISIDE SLIDER:
#  .heavisideSlider          Displays Heaviside and related functions
################################################################################


.heavisideSlider <- 
function()
{   
    # A function implemented by Diethelm Wuertz

    # Description
    #   Displays Heaviside and related functions

    # FUNCTION:
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        a = .sliderMenu(no = 1)
        by = .sliderMenu(no = 2)
        
        # Frame:
        par(mfrow = c(2, 2), cex = 0.7)
        
        # FGN TimeSeries:
        x = seq(-10, 10, by = by)
        H  = Heaviside(x, a)
        Sign  = Sign(x, a)
        Delta  = Delta(x, a)
        Boxcar  = Boxcar(x, a)
        Ramp  = Ramp(x, a)
        plot(x, H, type = "b", ylim = c(-1, 1), 
            col = "steelblue", pch = 19)
        title(main = paste("Heaviside | a =", a))
        grid()
        plot(x, Sign, type = "b", ylim = c(-1, 1), 
            col = "steelblue", pch = 19)
        title(main = paste("Sign | a =", a))
        grid()
        plot(x, Boxcar, type = "b", ylim = c(-1, 1), 
            col = "steelblue", pch = 19)
        title(main = paste("Boxcar | a =", a))
        grid()
        plot(x, Ramp, type = "b", 
            col = "steelblue", pch = 19)
        title(main = paste("Ramp | a =", a))
        grid()
        
        # Reset Frame:
        par(mfrow = c(1, 1), cex = 0.7)
    }
  
    # Open Slider Menu:
    .sliderMenu(refresh.code,
       names =       c(  "a", "by"),
       minima =      c(  -10,  0.1),
       maxima =      c(   10,  1.0),
       resolutions = c(  0.1,  0.1),
       starts =      c(    0,  0.2))
}


################################################################################

