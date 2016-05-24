#------------------------------------------------------------------------------#
#                     R Interface to C API of IBM ILOG CPLEX                   #
#------------------------------------------------------------------------------#

#  cplex_checkAPI.R
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
#                              the interface (check)                           #
#------------------------------------------------------------------------------#

checkCopyLpCPLEX <- function(env, lp, nCols, nRows, lpdir, objf, rhs, sense,
                             matbeg, matcnt, matind, matval, lb, ub,
                             rngval = NULL
                            ) {

    if (is.null(rngval)) {
        Crngval <- as.null(rngval)
    }
    else {
        Crngval <- as.numeric(rngval)
    }

    status <- .Call("checkCopyLp", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.integer(nCols),
                    as.integer(nRows),
                    as.integer(lpdir),
                    as.numeric(objf),
                    as.numeric(rhs),
                    as.character(paste(sense, collapse = "")),
                    as.integer(matbeg),
                    as.integer(matcnt),
                    as.integer(matind),
                    as.numeric(matval),
                    as.numeric(lb),
                    as.numeric(ub),
                    Crngval
              )

    return(status)
}


#------------------------------------------------------------------------------#

checkCopyLpwNamesCPLEX <- function(env, lp, nCols, nRows, lpdir, objf, rhs,
                                   sense, matbeg, matcnt, matind, matval,
                                   lb, ub, rngval = NULL,
                                   cnames = NULL, rnames = NULL) {

    if (is.null(rngval)) {
        Crngval <- as.null(rngval)
    }
    else {
        Crngval <- as.numeric(rngval)
    }

    if (is.null(cnames)) {
        Ccnames <- as.null(cnames)
    }
    else {
        Ccnames <- as.character(cnames)
    }

    if (is.null(rnames)) {
        Crnames <- as.null(rnames)
    }
    else {
        Crnames <- as.character(rnames)
    }

    status <- .Call("checkCopyLpwNames", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.integer(nCols),
                    as.integer(nRows),
                    as.integer(lpdir),
                    as.numeric(objf),
                    as.numeric(rhs),
                    as.character(paste(sense, collapse = "")),
                    as.integer(matbeg),
                    as.integer(matcnt),
                    as.integer(matind),
                    as.numeric(matval),
                    as.numeric(lb),
                    as.numeric(ub),
                    Crngval,
                    Ccnames,
                    Crnames
              )

    return(status)
}


#------------------------------------------------------------------------------#

checkCopyQuadCPLEX <- function(env, lp, qmatbeg, qmatcnt, qmatind, qmatval) {

    status <- .Call("checkCopyQuad", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.integer(qmatbeg),
                    as.integer(qmatcnt),
                    as.integer(qmatind),
                    as.numeric(qmatval)
              )

    return(status)
}


#------------------------------------------------------------------------------#

checkCopyQPsepCPLEX <- function(env, lp, qsepvec) {

    status <- .Call("checkCopyQPsep", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.numeric(qsepvec)
              )

    return(status)
}


#------------------------------------------------------------------------------#

checkAddRowsCPLEX <- function(env, lp, ncols, nrows, nnz,
                              matbeg, matind, matval,
                              rhs = NULL, sense = NULL,
                              cnames = NULL, rnames = NULL) {

    if (is.null(rhs)) {
        Crhs <- as.null(rhs)
    }
    else {
        Crhs <- as.numeric(rhs)
    }

    if (is.null(sense)) {
        Csense <- as.null(sense)
    }
    else {
        Csense <- as.character(paste(sense, collapse = ""))
    }

    if (is.null(cnames)) {
        Ccnames <- as.null(cnames)
    }
    else {
        Ccnames <- as.character(cnames)
    }

    if (is.null(rnames)) {
        Crnames <- as.null(rnames)
    }
    else {
        Crnames <- as.character(rnames)
    }

    status <- .Call("checkAddRows", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.integer(ncols),
                    as.integer(nrows),
                    as.integer(nnz),
                    Crhs,
                    Csense,
                    as.integer(matbeg),
                    as.integer(matind),
                    as.numeric(matval),
                    Ccnames,
                    Crnames
              )

    return(status)
}


#------------------------------------------------------------------------------#

checkAddColsCPLEX <- function(env, lp, ncols, nnz, objf, matbeg, matind, matval,
                              lb = NULL, ub = NULL, cnames = NULL) {

    if (is.null(lb)) {
        Clb <- as.null(lb)
    }
    else {
        Clb <- as.numeric(lb)
    }

    if (is.null(ub)) {
        Cub <- as.null(ub)
    }
    else {
        Cub <- as.numeric(ub)
    }

    if (is.null(cnames)) {
        Ccnames <- as.null(cnames)
    }
    else {
        Ccnames <- as.character(cnames)
    }

    status <- .Call("checkAddCols", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.integer(ncols),
                    as.integer(nnz),
                    as.numeric(objf),
                    as.integer(matbeg),
                    as.integer(matind),
                    as.numeric(matval),
                    Clb,
                    Cub,
                    Ccnames
              )

    return(status)
}


#------------------------------------------------------------------------------#

checkChgCoefListCPLEX <- function(env, lp, nnz, ia, ja, ra) {

    status <- .Call("checkChgCoefList", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.integer(nnz),
                    as.integer(ia),
                    as.integer(ja),
                    as.numeric(ra)
              )

    return(status)
}


#------------------------------------------------------------------------------#

checkCopyColTypeCPLEX <- function(env, lp, xctype) {

    status <- .Call("checkCopyColType", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.character(paste(xctype, collapse = ""))
              )

    return(status)
}


#------------------------------------------------------------------------------#

checkValsCPLEX <- function(env, lp, nval,
                           rind = NULL, cind = NULL, val = NULL) {

    if (is.null(rind)) {
        Crind <- as.null(rind)
    }
    else {
        Crind <- as.integer(rind)
    }

    if (is.null(cind)) {
        Ccind <- as.null(cind)
    }
    else {
        Ccind <- as.integer(cind)
    }

    if (is.null(val)) {
        Cval <- as.null(val)
    }
    else {
        Cval <- as.numeric(val)
    }

    status <- .Call("checkVals", PACKAGE = "cplexAPI",
                    cplexPointer(env),
                    cplexPointer(lp),
                    as.integer(nval),
                    Crind,
                    Ccind,
                    Cval
              )

    return(status)
}

