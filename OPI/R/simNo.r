#
# An implementation of the OPI that simulates a patient that never responds.
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: October 2012
#
# Copyright 2012 Andrew Turpin
# This program is part of the OPI (http://perimetry.org/OPI).
# OPI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

simNo.opiClose         <- function() { return(NULL) }
simNo.opiQueryDevice   <- function() { return (list(type="SimNo")) }

################################################################################
# Input
#   display Dimensions of plot area to display stim. c(-x,+x,-y,+y) No display if NULL
#
# Return NULL if succesful, string error message otherwise  
################################################################################
simNo.opiInitialize <- function(display=NULL) {
    if(simDisplay.setupDisplay(display))
        warning("opiInitialize (SimNo): display parameter may not contain 4 numbers.")
    return(NULL)
}

################################################################################
# Set background of plot area to col
# Return:
#   NULL - succsess
#   -1   - opiInitialize not called
################################################################################
simNo.opiSetBackground <- function(col, gridCol) { 
    return (simDisplay.setBackground(col, gridCol))
}

################################################################################
#
################################################################################
simNo.opiPresent <- function(stim) {
    simDisplay.present(stim$x, stim$y, stim$color, stim$duration, stim$responseWindow)

    return ( list(
        err = NULL,
        seen= FALSE,
        time= 0
    ))
}
