# Copyright 2007, 2008, 2010 Mario Pineda-Krch.
#
# This file is part of the R package GillespieSSA.
#
# GillespieSSA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# GillespieSSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GillespieSSA.  If not, see <http://www.gnu.org/licenses/>.

ssa.terminate <- function(args,out.rxn,tf,method,maxWallTime,verbose) {

  # Get the final time and state vector
  t <- out.rxn$timeSeries[dim(out.rxn$timeSeries)[1],1]
  x <- out.rxn$timeSeries[dim(out.rxn$timeSeries)[1],2:dim(out.rxn$timeSeries)[2]]

  # Figure out all the reasons why the simulation terminated
  terminationStatus <- NULL 
  if (t>=tf)          terminationStatus <- c(terminationStatus, "finalTime")
  if (all(x==0))      terminationStatus <- c(terminationStatus, "extinction")
  if (any(x<0))       terminationStatus <- c(terminationStatus, "negativeState")
  if (all(out.rxn$eval_a==0)) terminationStatus <- c(terminationStatus, "zeroProp")
  if (out.rxn$elapsedWallTime>=maxWallTime) terminationStatus <- c(terminationStatus, "maxWallTime") 

  # Calculate some stats for the used method
  stats <- list(startWallime       = out.rxn$startWallTime,
                endWallTime        = out.rxn$endWallTime,
                elapsedWallTime    = out.rxn$elapsedWallTime,
                terminationStatus  = terminationStatus,
                nSteps             = length(out.rxn$stepSize),
                meanStepSize       = mean(out.rxn$stepSize),
                sdStepSize         = sd(out.rxn$stepSize),
                nSuspendedTauLeaps = out.rxn$nSuspendedTauLeaps)

  # Figure out why the simulation terminated and print some info/stats
  if (verbose) {
    cat("tf: ",t,"\n","TerminationStatus: ",sep="")
    cat(stats$terminationStatus,sep=",")
    cat("\nDuration: ",stats$elapsedWallTime," seconds\n",
        "Method: ",method,"\n",
        "Nr of steps: ",stats$nSteps,"\n",
        "Mean step size: ",stats$meanStepSize,"+/-",stats$sdStepSize,"\n",sep="")
    if (method=="OTL") {
      cat("Nr suspended tau leaps: ",stats$nSuspendedTauLeaps,
          "(",100*(round(stats$nSuspendedTauLeaps/stats$nSteps)),"%)\n",sep="")
    }
    cat("End wall time: ",stats$endWallTime,"\n",sep="")
    cat("--------------------\n")
    }

  # Return simulation results ('chopping' off any rows in the timeSeries matrix that have no values (NA))    
  return(list(data  = out.rxn$timeSeries, 
              stats = stats, 
              args  = args))
}