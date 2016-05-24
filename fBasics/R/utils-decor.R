
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
# FUNCTION:        DESCRIPTION:
#  decor            Adds horizontal grid and L shaped box
#  hgrid            Adds horizontal grid lines
#  vgrid            Adds vertical grid lines
#  boxL             Adds L-shaped box
#  box_             Adds unterlined box
#  .xrug            Adds rugs on x axis
#  .yrug            Adds rugs on y axis
#  copyright        Adds copyright notice
################################################################################


decor <-
function()
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Adds a horizontal grid and L shaped box
   
    # FUNCTION:

    # Add:
    hgrid()
    boxL()
    
    # Return Value:
    invisible()
}


################################################################################


hgrid <-
function(ny = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds horizontal grid lines
    
    # FUNCTION:

    # Add:
    grid(NA, ny, ...)
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


vgrid <-
function(nx = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds vertical grid lines
    
    # FUNCTION:

    # Add:
    grid(nx, NA, ...)
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


boxL <-
function(col = "white")
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Adds L-shaped box
    
    # Add:
    box()
    box(bty = "7", col = col)
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


box_ <-
function(col = c("white", "black"))
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds unterlined box
    
    # FUNCTION:

    # Add:
    box(bty = "c", col = col[1])
    box(bty = "]", col = col[2])
    box(bty = "7", col = col[1])
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.xrug <-
function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds rugs on x axis
    
    # FUNCTION:

    # Add:
    rug(as.vector(x), ticksize = 0.01, side = 1, quiet = TRUE)
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.yrug <-
function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds rugs on y axis
    
    # FUNCTION:

    # Add:
    rug(as.vector(x), ticksize = 0.01, side = 2, quiet = TRUE)
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


copyright <-
function()
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds copyright notice
    
    # FUNCTION:

    # Add:
    Year = substr(Sys.Date(), 1, 4)
    mtext(paste("(c) Rmetrics", Year),
        side = 4, line = 0, adj = 0,
        font = 1, cex = 0.7*par("cex"), col = "grey")
        
    # Return Value:
    invisible()
}


################################################################################


