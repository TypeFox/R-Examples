##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

# This stub function simply calls down to a stub in the library.

ptest <- function()
{
    return_val <- .Call("ptest")

    # If the value is numeric then it means that
    # MPI is not initialized and the function should abort
    # and return FALSE
    if ( is.numeric(return_val) ) {
        warning(paste("MPI is not initialized. Function is aborted.\n"))
        return_val <- FALSE
    }

    return(return_val)
}
