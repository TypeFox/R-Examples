
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
# FUNCTION:             DESCRIPTION:
#  gedSlider             Displays Generalized Error Distribution and RVS
################################################################################


gedSlider <-  
function(type = c("dist", "rand"))
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Displays interactively skew GED distribution
    
    # Note:
    #   dged(x, mean = 0, sd = 1, nu = 5)
    
    # FUNCTION:
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        N      = .sliderMenu(no = 1)
        mean   = .sliderMenu(no = 2)
        sd     = .sliderMenu(no = 3)
        nu     = .sliderMenu(no = 4)
        
        # Compute Data:  
        xmin = round(qged(0.01, mean, sd, nu), digits = 2)
        xmax = round(qged(0.99, mean, sd, nu), digits = 2)
        s = seq(xmin, xmax, length = N)
        y1 = dged(s, mean, sd, nu)
        y2 = pged(s, mean, sd, nu)
        main1 = paste("GED Density\n", 
            "mean = ", as.character(mean), " | ",
            "sd = ", as.character(sd), " | ",
            "nu = ", as.character(nu) )
        main2 = paste("GED Probability\n",
            "xmin [0.01] = ", as.character(xmin), " | ",
            "xmax [0.99] = ", as.character(xmax) )   
            
        # Random Numbers:
        if (type[1] == "rand") {
            x = rged(N, mean, sd, nu)
        }      
        
        # Frame:    
        par(mfrow = c(2, 1), cex = 0.7)
        
        # Density:
        if (type[1] == "rand") {
            hist(x, probability = TRUE, col = "steelblue", border = "white",
                breaks = "FD",
                xlim = c(xmin, xmax), ylim = c(0, 1.1*max(y1)), main = main1 )
            lines(s, y1, col = "orange")
        } else {
            plot(s, y1, type = "l", xlim = c(xmin, xmax), col = "steelblue")
            abline (h = 0, lty = 3)
            title(main = main1)  
            grid()
        }
            
        # Probability:
        plot(s, y2, type = "l", xlim = c(xmin, xmax), ylim = c(0, 1),
            col = "steelblue" )
        abline (h = 0, lty = 3)
        title(main = main2) 
        grid()
        
        # Frame:
        par(mfrow = c(1, 1), cex = 0.7) 
    }
  
    # Open Slider Menu:
    .sliderMenu(refresh.code,
       names =       c(   "N", "mean",  "sd",  "nu"),
       minima =      c(   10,    -5.0,   0.1,   2.1),
       maxima =      c( 1000,    +5.0,   5.0,  10.0),
       resolutions = c(   10,     0.1,   0.1,   0.1),
       starts =      c(  100,     0.0,   1.0,   5.0)
    )
}


################################################################################

