
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
# FUNCTION:                    DESCRIPTION:
#   show.fPFOLIOBACKTEST        Print method for 'fPFOLIOBACKTEST' objects
################################################################################


setMethod("show", "fPFOLIOBACKTEST",
function(object)
{
    # A function implemented by Diethelm Wuertz and William Chen
    
    # Description:
    #   S4 Print Method for an object of class "fPFOLIODATA"

    # Arguments:
    #   object - an object of class "fPFOLIOSPEC"

    # FUNCTION:

    # Windows:
    cat("\nBacktest Specification:\t")
    cat("\n\n Windows Function:           ",
        object@windows$windows)
    cat("\n Windows Params:\t")
    winParams <- object@windows$params
    paramNames <- names(winParams)
    chars <- nchar(paramNames)
    for (i in seq(along = winParams)){
        cat("\n  -", paramNames[i], 
            paste(rep(" ", 24-nchar(paramNames[i])),collapse = ""), 
            as.character(winParams[[i]]))
    }
 
    # Strategy:
    cat("\n\n Strategy Function:          ",
        object@strategy$strategy)
    cat("\n Strategy Params:\t")
    strategyParams <- object@strategy$params
    paramNames <- names(strategyParams)
    for (i in seq(along = strategyParams)){
        cat("\n  -", paramNames[i], 
        paste(rep(" ", 24-nchar(paramNames[i])),collapse = ""), 
        as.character(strategyParams[[i]]))      
    }     
    
    # Smoother:
    cat("\n\n Smoother Function:          ",
        object@smoother$smoother)
    cat("\n Smoother Params:\t")
    smootherParams <- object@smoother$params
    paramNames <- names(smootherParams)
    for (i in seq(along = smootherParams)){
        cat("\n  -", paramNames[i], 
        paste(rep(" ", 24-nchar(paramNames[i])),collapse = ""), 
        substr(as.character(smootherParams[[i]]), 1, 5))      
    }     
    
    # Messages:
    cat("\n\n Messages:                   ",
        unlist(object@messages), "\n\n")
   
    # Return Value:
    invisible(object)
})


################################################################################

