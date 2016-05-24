#------------------------------------------------------------------------------#
#                           R interface to COIN-OR Clp                         #
#------------------------------------------------------------------------------#

#  clp.R
#  R interface to COIN-OR Clp.
#
#  Copyright (C) 2011-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of clpAPI.
#
#  ClpAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ClpAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with clpAPI  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#

# return codes of clp optimizations
return_codeCLP <- function(code) {
    if (code == 0) { return( "solution process was successful" ) }
    else { return(paste("Failed to obtain solution, unknown error code:", code)) }
}


# solution status codes of clp optimizations
status_codeCLP <- function(code) {
    if (code == 0)      { return( "solution is optimal" ) }
    else if (code == 1) { return( "solution is primal infeasible" ) }
    else if (code == 2) { return( "solution is dual infeasible" ) }
    else if (code == 3) { return( "stopped on iterations etc" ) }
    else if (code == 4) { return( "stopped due to errors" ) }
    else { return(paste("unknown status code:", code)) }
}


