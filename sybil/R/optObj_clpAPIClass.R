#  optObj_clpAPIClass.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#  
#  This file is part of sybil.
#
#  Sybil is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Sybil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybil.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                   definition of the class optObj_clpAPI                      #
#------------------------------------------------------------------------------#

setClass(Class = "optObj_clpAPI", contains = "optObj")


#------------------------------------------------------------------------------#
#                                  methods                                     #
#------------------------------------------------------------------------------#

setMethod("delProb", signature(lp = "optObj_clpAPI"),

    function(lp, ...) {

        clpAPI::delProbCLP(lp@oobj)

    }
)


#------------------------------------------------------------------------------#

setMethod("initProb", signature(lp = "optObj_clpAPI"),

    function(lp, to = NULL, ...) {

        lp@oobj <- clpAPI::initProbCLP()

        if (is.null(to)) {
            clpAPI::setLogLevelCLP(lp@oobj, 0)
        }
        else {
            stopifnot(is(to, "numeric"))
            clpAPI::setLogLevelCLP(lp@oobj, to)
        }

        return(lp)
    }
)


#------------------------------------------------------------------------------#

setMethod("backupProb", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- FALSE
        np  <- FALSE

        fname <- tempfile(pattern = "CLP_PROB_", fileext = ".tmp")
        ft <- clpAPI::saveModelCLP(lp@oobj, fname)
        if (ft != 0) {
            stop("cannot save model")
        }

        # create new lp problem
        np <- clpAPI::initProbCLP()
        clpAPI::setLogLevelCLP(np, 0)
        ft <- clpAPI::restoreModelCLP(np, fname)
        if (ft != 0) {
            stop("cannot read model")
        }
        else {
            unlink(fname)
        }

        # create new optObj object
        if (!identical(np, FALSE)) {
            out <- new("optObj_clpAPI", lp@solver, lp@method, lp@probType)
            out@oobj <- np
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#


setMethod("setSolverParm", signature(lp = "optObj_clpAPI"),

    function(lp, solverParm) {

        out <- FALSE

        wrong_solver_msg(lp, "setSolverParm")

#        if ( ! ((is.data.frame(solverParm)) || (is.list(solverParm))) ) {
#            stop(sQuote(solverParm), " must be list or data.frame")
#        }
#
#        if (any(is.na(solverParm))) {
#            stop(sQuote(solverParm), " contains NA values")
#        }

        # no parameters in COIN-OR CLP yet.
        #    lp@oobj <- clpAPI::initProbCLP()
        #    clpAPI::setLogLevelCLP(lp@oobj, 0)

        return(out)

    }
)


#------------------------------------------------------------------------------#

setMethod("getSolverParm", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- FALSE

        wrong_solver_msg(lp, "getSolverParm")

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_clpAPI", lpdir = "character"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == "max", -1, 1)
        clpAPI::setObjDirCLP(lp@oobj, dr)

    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_clpAPI", lpdir = "numeric"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == -1, -1, 1)
        clpAPI::setObjDirCLP(lp@oobj, dr)

    }
)

#------------------------------------------------------------------------------#

setMethod("getObjDir", signature(lp = "optObj_clpAPI"),

    function(lp) {

        dr <- clpAPI::getObjDirCLP(lp@oobj)
        if (dr == -1) {
            out <- "max"
        }
        else if (dr == 1) {
            out <- "min"
        }
        else {
            out <- FALSE
        }
        return(out)

    }
)


#------------------------------------------------------------------------------#

setMethod("addRows", signature(lp = "optObj_clpAPI", nrows = "numeric"),

    function(lp, nrows) {

        ncols <- clpAPI::getNumColsCLP(lp@oobj)
        out   <- clpAPI::resizeCLP(lp@oobj, nrows, ncols)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addCols", signature(lp = "optObj_clpAPI", ncols = "numeric"),

    function(lp, ncols) {

        nrows <- clpAPI::getNumRowsCLP(lp@oobj)
        out   <- clpAPI::resizeCLP(lp@oobj, nrows, ncols)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsCols", signature(lp = "optObj_clpAPI",
                                   nrows = "numeric", ncols = "numeric"),

    function(lp, nrows, ncols) {

        out <- clpAPI::resizeCLP(lp@oobj, nrows, ncols)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumRows", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- clpAPI::getNumRowsCLP(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumCols", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- clpAPI::getNumColsCLP(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addColsToProb", signature(lp = "optObj_clpAPI"),

    # j: vector containing the new column indices (must be ascending)
    # rind: list, containing the row indices of the new nz elements
    # nzval: list, containing the new nz elements
    #
    # j, obj, lb, rind and nzval must have the same length

    function(lp, j, obj, lb, ub, rind, nzval) {

        cst <- c(0, cumsum(unlist(lapply(rind, length))))
        print(cst)
        print(unlist(rind)-1)
        print(unlist(nzval))
        out <- clpAPI::addColsCLP(lp@oobj, length(j), lb, ub, obj,
                                  cst, unlist(rind)-1, unlist(nzval))

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsToProb", signature(lp = "optObj_clpAPI"),

    # i: vector containing the new row indices (must be ascending)
    # cind: list, containing the column indices of the new nz elements
    # nzval: list, containing the new nz elements
    #
    # i, type, lb, ub, cind and nzval must have the same length
    #
    # type can be one of the following:
    # "F" = free variable                -INF <  x <  INF
    # "L" = variable with lower bound      lb <= x <  INF
    # "U" = variable with upper bound    -INF <  x <= ub
    # "D" = double-bounded variable        lb <= x <= ub
    # "E" = fixed variable                 lb  = x  = ub

    function(lp, i, type, lb, ub, cind, nzval, rnames = NULL) {

        stopifnot(length(lb) == length(ub))
        cub <- ub
        ebc <- type == "E"
        cub[ebc] <- lb[ebc]

        cst <- c(0, cumsum(unlist(lapply(cind, length))))

        out <- clpAPI::addRowsCLP(lp@oobj, length(i), lb, cub,
                                  cst, unlist(cind)-1, unlist(nzval))

        #if (!is.null(rnames)) {
        #    for (rind in seq(along = i)) {
        #        clpAPI::rowNameCLP(lp@oobj, i = i[rind], rname = rnames[rind])
        #    }
        #}

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBnds", signature(lp = "optObj_clpAPI"),

    function(lp, j, lb, ub) {

        tmp_lb <- clpAPI::getColLowerCLP(lp@oobj)
        tmp_ub <- clpAPI::getColUpperCLP(lp@oobj)
        tmp_lb[j] <- lb
        tmp_ub[j] <- ub
        clpAPI::chgColLowerCLP(lp@oobj, tmp_lb)
        clpAPI::chgColUpperCLP(lp@oobj, tmp_ub)

    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBndsObjCoefs", signature(lp = "optObj_clpAPI"),

    function(lp, j, lb, ub, obj_coef) {

        # usable only for model creation!
        clpAPI::chgColLowerCLP(lp@oobj, lb)
        clpAPI::chgColUpperCLP(lp@oobj, ub)
        clpAPI::chgObjCoefsCLP(lp@oobj, obj_coef)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsLowBnds", signature(lp = "optObj_clpAPI", j = "numeric"),

    function(lp, j) {

        out <- clpAPI::getColLowerCLP(lp@oobj)[j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsUppBnds", signature(lp = "optObj_clpAPI", j = "numeric"),

    function(lp, j) {

        out <- clpAPI::getColUpperCLP(lp@oobj)[j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeRowsBnds", signature(lp = "optObj_clpAPI"),

    function(lp, i, lb, ub) {

        tmp_lb <- clpAPI::getRowLowerCLP(lp@oobj)
        tmp_ub <- clpAPI::getRowUpperCLP(lp@oobj)
        tmp_lb[i] <- lb
        tmp_ub[i] <- ub
        clpAPI::chgRowLowerCLP(lp@oobj, tmp_lb)
        clpAPI::chgRowUpperCLP(lp@oobj, tmp_ub)
    }
)


#------------------------------------------------------------------------------#

setMethod("setRhsZero", signature(lp = "optObj_clpAPI"),

    function(lp) {

        nrows <- clpAPI::getNumRowsCLP(lp@oobj)
        zeros <- rep(0, nrows)
        clpAPI::chgRowLowerCLP(lp@oobj, zeros)
        clpAPI::chgRowUpperCLP(lp@oobj, zeros)

    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsLowBnds", signature(lp = "optObj_clpAPI", i = "numeric"),

    function(lp, i) {

        out <- clpAPI::getRowLowerCLP(lp@oobj)[i]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsUppBnds", signature(lp = "optObj_clpAPI", i = "numeric"),

    function(lp, i) {

        out <- clpAPI::getRowUpperCLP(lp@oobj)[i]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeObjCoefs", signature(lp = "optObj_clpAPI"),

    function(lp, j, obj_coef) {

        tmp_obj_coef <- clpAPI::getObjCoefsCLP(lp@oobj)
        tmp_obj_coef[j] <- obj_coef
        clpAPI::chgObjCoefsCLP(lp@oobj, tmp_obj_coef)

    }
)


#------------------------------------------------------------------------------#

setMethod("getObjCoefs", signature(lp = "optObj_clpAPI", j = "numeric"),

    function(lp, j) {

        out <- clpAPI::getObjCoefsCLP(lp@oobj)[j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("loadLPprob", signature(lp = "optObj_clpAPI"),

    function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype,
             lpdir = "max", rub = NULL, ctype = NULL,
             cnames = NULL, rnames = NULL, pname = NULL,
             defLowerBnd = SYBIL_SETTINGS("MAXIMUM") * -1,
             defUpperBnd = SYBIL_SETTINGS("MAXIMUM")) {

        stopifnot(is(mat, "Matrix"))
        
        if (!is.null(ctype)) {
            warning("argument ", sQuote(ctype), " is currently ignored")
        }

        if (is.null(rub)) {
            stopifnot(is(defLowerBnd, "numeric"), is(defUpperBnd, "numeric"))
            ##crub <- numeric(nRows)
            # Default value for upper bound is SYBIL_SETTINGS("MAXIMUM").
            
            # If a constraint has an upper bound (type "U"), the corresponding
            # value is copied from rlb to crub and crlb get's the default
            # lower bound SYBIL_SETTINGS("MAXIMUM") * -1.
            
            # If a constraint is equality (type "E"), the corresponding
            # value is copied from rlb to crub.

            crub <- rep(defUpperBnd, nRows)
            crlb <- rlb
            ubc  <- rtype == "U"
            crub[ubc] <- rlb[ubc]
            crlb[ubc] <- defLowerBnd
            ebc  <- rtype == "E"
            crub[ebc] <- rlb[ebc]
        }
        else {
            stopifnot(length(rub) == length(rlb))
            crub <- rub
            crlb <- rlb
            ebc  <- rtype == "E"
            crub[ebc] <- crlb[ebc]
        }

        # optimization direction
        setObjDir(lp, lpdir = lpdir)
        
        # load problem
        TMPmat <- as(mat, "CsparseMatrix")
        clpAPI::loadProblemCLP(lp@oobj,
                               ncols    = nCols,
                               nrows    = nRows,
                               ia       = TMPmat@i,
                               ja       = TMPmat@p,
                               ra       = TMPmat@x,
                               lb       = lb,
                               ub       = ub,
                               obj_coef = obj,
                               rlb      = crlb,
                               rub      = crub)

        # row and column names
        if ( (!is.null(rnames)) && (!is.null(cnames)) ) {
            clpAPI::copyNamesCLP(lp@oobj, cnames = cnames, rnames = rnames)
        }

        # problem name
        if (!is.null(rnames)) {
            clpAPI::probNameCLP(lp@oobj, pname = pname)
        }
    }
)


#------------------------------------------------------------------------------#

setMethod("scaleProb", signature(lp = "optObj_clpAPI"),

    function(lp, opt) {

        clpAPI::scaleModelCLP(lp@oobj, opt)

    }
)


#------------------------------------------------------------------------------#

setMethod("solveLp", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- FALSE
        switch(lp@method,
            "inidual" = {
                out <- clpAPI::solveInitialDualCLP(lp@oobj)
            },
            "iniprimal" = {
                out <- clpAPI::solveInitialPrimalCLP(lp@oobj)
            },
            "inibarrier" = {
                out <- clpAPI::solveInitialBarrierCLP(lp@oobj)
            },
            "inibarriernoc" = {
                out <- clpAPI::solveInitialBarrierNoCrossCLP(lp@oobj)
            },
            "dual" = {
                out <- clpAPI::dualCLP(lp@oobj)
            },
            "primal" = {
                out <- clpAPI::primalCLP(lp@oobj)
            },
            "idiot" = {
                clpAPI::idiotCLP(lp@oobj) # idiotCLP has no return value
                out <- 0
            },
            {
                out <- clpAPI::solveInitialCLP(lp@oobj)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjVal", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- clpAPI::getObjValCLP(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRedCosts", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- clpAPI::getColDualCLP(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolStat", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- clpAPI::getSolStatusCLP(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getFluxDist", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- clpAPI::getColPrimCLP(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColPrim", signature(lp = "optObj_clpAPI", j = "numeric"),

    function(lp, j) {

        out <- clpAPI::getColPrimCLP(lp@oobj)[j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumNnz", signature(lp = "optObj_clpAPI"),

    function(lp) {

        out <- clpAPI::getNumNnzCLP(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("writeProb", signature(lp = "optObj_clpAPI", fname = "character"),

    function(lp, fname, ff = "lp") {

        out <- clpAPI::saveModelCLP(lp@oobj, fname = fname)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("readProb", signature(lp = "optObj_clpAPI", fname = "character"),

    function(lp, fname, ff = "mps", ...) {

        switch(ff,
            "mps" = {
                fl <- clpAPI::readMPSCLP(lp@oobj, fname = fname, ...)
            },
            "clp" = {
                fl <- clpAPI::restoreModelCLP(lp@oobj, fname = fname)
            },
            {
                message("wrong format!")
                fl <- 1
            }
        )
        out <- ifelse(fl == 0, TRUE, fl)
        
        return(out)
    }
)


#------------------------------------------------------------------------------#

