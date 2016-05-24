#------------------------------------------------------------------------------#
#                     R Interface to C API of IBM ILOG CPLEX                   #
#------------------------------------------------------------------------------#

#  cplex.R
#  R Interface to C API of IBM ILOG CPLEX Version 12.1 to 12.6.
#
#  Copyright (C) 2011-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of cplexAPI.
#
#  CplexAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  CplexAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with cplexAPI.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
# for checkOptSol

# return codes of ILOG CPLEX optimizations
return_codeCPLEX <- function(code) {

    if (code == 0) {
        error <- "Optimization ended successfull"
    } else {
        error <- getErrorStrCPLEX(code)
    }

    return(error)

}


# status codes of ILOG CPLEX optimizations
status_codeCPLEX <- function(env, code) {

    return(getStatStrCPLEX(env, code))

}


#------------------------------------------------------------------------------#
# get values and names of non-default parameters
getParmValCPLEX <- function(env) {

    parmCname <- getChgParmCPLEX(env)

    if (!is.null(parmCname)) {

        parmType     <- sapply(parmCname, function(x) getParmTypeCPLEX(env, x))
        parmCnameStr <- sapply(parmCname, function(x) getParmNameCPLEX(env, x))

        intVal <- sapply(parmCname[parmType == CPX_PARAMTYPE_INT],
                         function(x) getIntParmCPLEX(env, x)
                        )
        dblVal <- sapply(parmCname[parmType == CPX_PARAMTYPE_DOUBLE],
                         function(x) getDblParmCPLEX(env, x)
                        )
        strVal <- sapply(parmCname[parmType == CPX_PARAMTYPE_STRING],
                         function(x) getStrParmCPLEX(env, x)
                        )

        names(intVal) <- parmCnameStr[parmType == CPX_PARAMTYPE_INT]
        names(dblVal) <- parmCnameStr[parmType == CPX_PARAMTYPE_DOUBLE]
        names(strVal) <- parmCnameStr[parmType == CPX_PARAMTYPE_STRING]

        parms <- list(integer = intVal, double = dblVal, string = strVal)

    }
    else {
        parms <- parmCname
    }

    return(parms)

}



