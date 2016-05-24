
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
# FUNCTION:              DESCRIPTION:
#  snormSlider            Displays Normal Distribution and RVS
################################################################################


snormSlider <- 
function(type = c("dist", "rand"))
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Displays interactively skew Normal distribution
    
    # Note:
    #   dsnorm(x, mean = 0, sd = 1, xi = 1.5)
    
    # FUNCTION:
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        N      = .sliderMenu(no = 1)
        mean   = .sliderMenu(no = 2)
        sd     = .sliderMenu(no = 3)
        xi     = .sliderMenu(no = 4)
        invert = .sliderMenu(no = 5)
        
        # Compute Data:  
        if (invert == 1) xi = 1/xi
        xmin = round(qsnorm(0.001, mean, sd, xi), digits = 2)
        xmax = round(qsnorm(0.999, mean, sd, xi), digits = 2)
        s = seq(xmin, xmax, length = N)
        y1 = dsnorm(s, mean, sd, xi)
        y2 = psnorm(s, mean, sd, xi)
        main1 = paste("Skew Normal Density\n", 
            "mean = ", as.character(mean), " | ",
            "sd = ", as.character(sd), " | ",
            "xi = ", as.character(xi) )
        main2 = paste("Skew Normal Probability\n",
            "xmin [0.001] = ", as.character(xmin), " | ",
            "xmax [0.999] = ", as.character(xmax) ) 
            
        # Random Numbers:
        if (type[1] == "rand") {
            x = rsnorm(N, mean, sd, xi) 
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
       names =       c(   "N", "mean",   "sd",  "xi", "xi.inv"),
       minima =      c(   10,    -5.0,    0.1,   1.0,       0 ),
       maxima =      c(  500,    +5.0,    5.0,  10.0,       1 ),
       resolutions = c(   10,     0.1,    0.1,   0.1,       1 ),
       starts =      c(  100,     0.0,    1.0,   1.0,       0 )
    )
}


################################################################################

