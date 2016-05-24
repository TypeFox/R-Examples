
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################
# FUNCTION:                     DESCRIPTION:
#  .amplExec                     Executes AMPL run file for a given project
#  .amplExample                  Optimizes mean variance portfolio example
###############################################################################


.amplExec <- 
  function(project="ampl")
  {
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #    Executes AMPL run file for a given project.
    
    # Details:
    #    Note the following files must exist:
    #    [project].mod - the AMPL model file
    #    [project].dat - the AMPL data file
    #    [project].run - the AMPL run file
    
    # FUNCTION:
    
    # Execute run file:
    command <- paste("ampl -t -vs", paste(project, "run", sep="."))
    solve <- system(command, intern=TRUE)
    
    # Print:
    cat(solve, sep="\n")
    
    # Return Value:
    invisible(solve)
  }


# -----------------------------------------------------------------------------

.amplExample <- 
  function()
  {
    # A function Implemented by Diethelm Wuertz
    
    # Description:
    #    Optimizes mean variance portfolio example.
    
    # FUNCTION:
    
    # Load Dataset
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    
    # Portfolio Data:
    nAssets <- 6
    data <- 100 * LPP2005REC[, 1:nAssets]
    
    # Optimization Arguments: 
    objective <- list(dvec=rep(0, nAssets), Dmat=cov(data))
    lower <- 0
    upper <- 1
    linCons <- list(
      mat = rbind(
        budget = rep(1, times=nAssets), 
        returns = colMeans(data)),
      lower = c(
        budget = 1, 
        return = mean(data)),
      upper = c(
        budget = 1, 
        return = mean(data)))
    control <- list()
    
    # Default - AMPL Interface:
    ampl <- ramplQP(objective, lower, upper, linCons, control)
    
    # Return Value:
    ampl
  }


###############################################################################

