
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
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                DESCRIPTION:
#  .counts1Plot              Creates a plot of counts
#  .counts2Plot              Creates an alternative plot of counts
################################################################################

  
.counts1Plot= 
function(x, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a plot of counts
    
    # Arguments:
    #   x - a vector of counts
    
    # Note:
    #   From R-help, ref. Gabor Grothendiek
    
    # Example:
    #   x <- c(26,26,27,27,27,27,28,28,28,28,28,28,28,28,28,28,28,28,
    #       28,29,29,29,29,29,29,29,29,29,29,29,30,30,30,30,30,30,30,
    #       30,31,31,31,31,32,33,33,33,34,34,35,43); .counts1Plot(x)
    
    # FUNCTION:
    
    # Plot Range:
    xmin = min(x)
    xmax = max(x)
    xlim = c(xmin - 0.1*(xmax-xmin), xmax + 0.1*(xmax-xmin))

    # Plot
    plot(seq(y) - match(y, y) ~ y, list(y = sort(x)), xlim = xlim, 
        ylab = "", yaxt = "n", pch = 19, frame.plot = FALSE, ...)
        
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.counts2Plot <- 
function(x, ...) 
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a plot of counts
    
    # Arguments:
    #   x - a vector of counts
    
    # Note:
    #   From R-help, ref. Jinsong Zhao
    
    # Example:
    #   x <- c(26,26,27,27,27,27,28,28,28,28,28,28,28,28,28,28,28,28,
    #       28,29,29,29,29,29,29,29,29,29,29,29,30,30,30,30,30,30,30,
    #       30,31,31,31,31,32,33,33,33,34,34,35,43); .counts2Plot(x)
    
    # FUNCTION:
    
    # Tabulate:
    x.table <- table(x)
    x.x <- as.numeric(names(x.table))
    x.y <- as.numeric(x.table)
    n <- length(x.x)
    x.final <- NULL
    for (i in 1:n) {
        tmp <- data.frame(x = rep(x.x[i], x.y[i]), y = 1:x.y[i])
        x.final <- rbind(x.final, tmp)
    }
    
    # Plot Range:
    xmin = min(x)
    xmax = max(x)
    xlim = c(xmin - 0.1*(xmax-xmin), xmax + 0.1*(xmax-xmin))
    
    # Plot:
    plot(x.final, yaxt = "n", ylab = "", xlim = xlim, pch = 19, 
        frame.plot = FALSE, ...)
        
    # Return Value:
    invisible(x)
}

 
################################################################################

