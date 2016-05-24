#  optObj_cplexAPIClass.R
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
#                  definition of the class optObj_cplexAPI                     #
#------------------------------------------------------------------------------#

setClass(Class = "optObj_cplexAPI", contains = "optObj")


#------------------------------------------------------------------------------#
#                                  methods                                     #
#------------------------------------------------------------------------------#

setMethod("delProb", signature(lp = "optObj_cplexAPI"),

    function(lp, closeEnv = TRUE) {

        if (isTRUE(closeEnv)) {
            cplexAPI::closeProbCPLEX(list(env = lp@oobj@env,
                                          lp = lp@oobj@lp))
        } else {
            cplexAPI::delProbCPLEX(lp@oobj@env, lp@oobj@lp)
        }

    }
)


#------------------------------------------------------------------------------#

setMethod("initProb", signature(lp = "optObj_cplexAPI"),

    function(lp, to = FALSE, ...) {

        #lp@oobj <- cplexAPI::openProbCPLEX()
        tmp <- cplexAPI::openProbCPLEX()
        lp@oobj <- new("cplexPointer",
                       en = tmp[["env"]],
                       pr = tmp[["lp"]])

        if (is.null(to)) {
            too <- FALSE
        }
        else {
            too <- to
        }
        
        if (isTRUE(too)) {
            cplexAPI::setIntParmCPLEX(lp@oobj@env,
                                      cplexAPI::CPX_PARAM_SCRIND,
                                      cplexAPI::CPX_ON)
        }
        else {
            cplexAPI::setIntParmCPLEX(lp@oobj@env,
                                      cplexAPI::CPX_PARAM_SCRIND,
                                      cplexAPI::CPX_OFF)
        }

        return(lp)
    }
)


#------------------------------------------------------------------------------#

setMethod("backupProb", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        out <- FALSE
        np  <- FALSE

        np <- cplexAPI::cloneProbCPLEX(lp@oobj@env, lp@oobj@lp)

        # create new optObj object
        if (!identical(np, FALSE)) {
            out <- new("optObj_cplexAPI", lp@solver, lp@method, lp@probType)
            out@oobj <- new("cplexPointer", en = lp@oobj@env, pr = np)
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#


setMethod("setSolverParm", signature(lp = "optObj_cplexAPI"),

    function(lp, solverParm) {

        out <- FALSE

        if ( ! ((is.data.frame(solverParm)) || (is.list(solverParm))) ) {
            stop(sQuote(solverParm), " must be list or data.frame")
        }

        if (any(is.na(solverParm))) {
            stop(sQuote(solverParm), " contains NA values")
        }

        intdbl  <- sapply(solverParm, is.integer)
        strparm <- sapply(solverParm, is.numeric)
        int  <- solverParm[intdbl]
        dbl  <- solverParm[intdbl == FALSE & strparm == TRUE]
        char <- solverParm[strparm == FALSE]

        if (length(int) > 0) {
            intp <- sapply(names(int), function(x) eval(parse(text = x)))
            intv <- unlist(int)
            for (i in seq(along = int)) {
                out  <- cplexAPI::setIntParmCPLEX(lp@oobj@env, intp[i], intv[i])
            }
        }

        if (length(dbl) > 0) {
            dblp <- sapply(names(dbl), function(x) eval(parse(text = x)))
            dblv <- unlist(dbl)
            for (i in seq(along = dbl)) {
                out  <- cplexAPI::setDblParmCPLEX(lp@oobj@env, dblp[i], dblv[i])
            }
        }

        if (length(char) > 0) {
            charp <- sapply(names(char), function(x) eval(parse(text = x)))
            charv <- unlist(char)
            for (i in seq(along = char)) {
                out  <- cplexAPI::setStrParmCPLEX(lp@oobj@env,
                                                  charp[i], charv[i])
            }
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolverParm", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        out <- cplexAPI::writeParmCPLEX(lp@oobj@env,
                                        "cplex_parameters.prm")
        message("Wrote the file 'cplex_parameters.prm'.")

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_cplexAPI", lpdir = "character"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == "max",
                     cplexAPI::CPX_MAX,
                     cplexAPI::CPX_MIN)
        cplexAPI::setObjDirCPLEX(lp@oobj@env, lp@oobj@lp, dr)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_cplexAPI", lpdir = "integer"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == cplexAPI::CPX_MAX,
                     cplexAPI::CPX_MAX,
                     cplexAPI::CPX_MIN)
        cplexAPI::setObjDirCPLEX(lp@oobj@env, lp@oobj@lp, dr)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_cplexAPI", lpdir = "numeric"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == -1,
                     cplexAPI::CPX_MAX,
                     cplexAPI::CPX_MIN)
        cplexAPI::setObjDirCPLEX(lp@oobj@env, lp@oobj@lp, dr)
    }
)

#------------------------------------------------------------------------------#

setMethod("getObjDir", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        dr <- cplexAPI::getObjDirCPLEX(lp@oobj@env, lp@oobj@lp)
        if (dr == cplexAPI::CPX_MAX) {
            out <- "max"
        }
        else if (dr == cplexAPI::CPX_MIN) {
            out <- "min"
        }
        else {
            out <- FALSE
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRows", signature(lp = "optObj_cplexAPI", nrows = "numeric"),

    function(lp, nrows) {

        out <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp, nrows)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addCols", signature(lp = "optObj_cplexAPI", ncols = "numeric"),

    function(lp, ncols) {

        out <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp, ncols)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsCols", signature(lp = "optObj_cplexAPI",
                                   nrows = "numeric", ncols = "numeric"),

    function(lp, nrows, ncols) {

        outi <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp, nrows)
        outj <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp, ncols)
        out  <- c(outi, outj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumRows", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        out <- cplexAPI::getNumRowsCPLEX(lp@oobj@env, lp@oobj@lp)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumCols", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        out <- cplexAPI::getNumColsCPLEX(lp@oobj@env, lp@oobj@lp)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addColsToProb", signature(lp = "optObj_cplexAPI"),

    # j: vector containing the new column indices (must be ascending)
    # rind: list, containing the row indices of the new nz elements
    # nzval: list, containing the new nz elements
    #
    # j, obj, lb, rind and nzval must have the same length

    function(lp, j, obj, lb, ub, rind, nzval) {

        beg <- c(0, cumsum(unlist(lapply(rind, length))))
        out <- cplexAPI::addColsCPLEX(lp@oobj@env, lp@oobj@lp,
                                      length(j), length(unlist(nzval)), obj,
                                      beg, unlist(rind)-1, unlist(nzval),
                                      lb, ub)

        return(out)
    }
)


#------------------------------------------------------------------------------#

#setMethod("addRowsToProb", signature(lp = "optObj_cplexAPI"),
#
#    # i: vector containing the new row indices (must be ascending)
#    # cind: list, containing the column indices of the new nz elements
#    # nzval: list, containing the new nz elements
#    #
#    # i, type, lb, cind and nzval must have the same length
#    #
#    # type can be one of the following:
#    # "F" = free variable                -INF <  x <  INF
#    # "L" = variable with lower bound      lb <= x <  INF
#    # "U" = variable with upper bound    -INF <  x <= ub
#    # "D" = double-bounded variable        lb <= x <= ub
#    # "E" = fixed variable                 lb  = x  = ub
#    # "R" = ranged constraint
#
#    function(lp, i, type, lb, ub, cind, nzval, rnames = NULL) {
#
#        cptype = character(length(type))
#        for (l in seq(along = type)) {
#            cptype[l] <- switch(type[l],
#                "L" = { "G" },
#                "U" = { "L" },
#                "E" = { "E" },
#                "R" = { "R" },
#                      { "E" }
#            )
#        }
#
#        beg <- c(0, cumsum(unlist(lapply(cind, length))))
#        out <- cplexAPI::addRowsCPLEX(env = lp@oobj@env, lp = lp@oobj@lp,
#                                      ncols = 0, nrows = length(i),
#                                      nnz = length(unlist(nzval)),
#                                      matbeg = beg, matind = unlist(cind)-1,
#                                      matval = unlist(nzval), rhs = lb,
#                                      sense = cptype, rnames = rnames)
#
#        return(out)
#    }
#)


#------------------------------------------------------------------------------#

setMethod("addRowsToProb", signature(lp = "optObj_cplexAPI"),

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

        cptype = character(length(type))
        for (l in seq(along = type)) {
            cptype[l] <- switch(EXPR = type[l],
                "L" = { "G" },
                "U" = { "L" },
                "D" = { "R" },
                "E" = { "E" },
                      { "E" }
            )
        }

        stopifnot(length(lb) == length(ub))
        rng      <- cptype == "R"
        cub      <- abs(ub[rng] - lb[rng])     # range
        
        ubc      <- cptype == "L"
        clb      <- lb
        clb[ubc] <- ub[ubc]
        
        beg <- c(0, cumsum(unlist(lapply(cind, length))))
        out <- cplexAPI::addRowsCPLEX(env = lp@oobj@env, lp = lp@oobj@lp,
                                      ncols = 0, nrows = length(i),
                                      nnz = length(unlist(nzval)),
                                      matbeg = beg, matind = unlist(cind)-1,
                                      matval = unlist(nzval), rhs = clb,
                                      sense = cptype, rnames = rnames)

        # set ranged (double bounded constraints)
        if (sum(rng) > 0) {
            #print(i[rng])
            cplexAPI::chgRngValCPLEX(env = lp@oobj@env, lp = lp@oobj@lp,
                                     nrows = sum(rng),
                                     ind = i[rng]-1,
                                     val = cub)
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBnds", signature(lp = "optObj_cplexAPI"),

    function(lp, j, lb, ub) {

        out <- cplexAPI::chgColsBndsCPLEX(lp@oobj@env,
                                          lp@oobj@lp, j-1, lb, ub)
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBndsObjCoefs", signature(lp = "optObj_cplexAPI"),

    function(lp, j, lb, ub, obj_coef) {

        outb <- cplexAPI::chgColsBndsCPLEX(lp@oobj@env,
                                           lp@oobj@lp, j-1, lb, ub)
        outo <- cplexAPI::chgObjCPLEX(lp@oobj@env, lp@oobj@lp,
                                      length(j), j-1, obj_coef)
        out  <- c(outb, outo)
        # usable only for model creation!
        # out <- cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp,
        #                               length(j), obj_coef, lb, ub)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsLowBnds", signature(lp = "optObj_cplexAPI", j = "numeric"),

    function(lp, j) {

        out <- cplexAPI::getLowBndsIdsCPLEX(lp@oobj@env, lp@oobj@lp, j-1)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsUppBnds", signature(lp = "optObj_cplexAPI", j = "numeric"),

    function(lp, j) {

        out <- cplexAPI::getUppBndsIdsCPLEX(lp@oobj@env, lp@oobj@lp, j-1)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeRowsBnds", signature(lp = "optObj_cplexAPI"),

    function(lp, i, lb, ub) {

#        out <- cplexAPI::chgRhsCPLEX(lp@oobj@env, lp@oobj@lp,
#                                     length(i), i-1, lb)

        stopifnot(length(lb) == length(ub))
        clb <- lb
       
        ct <- mapply(cplexAPI::getSenseCPLEX, i-1, i-1,
                     MoreArgs = list(env = lp@oobj@env, lp  = lp@oobj@lp))
        
        # If a constraint is a ranged constraint, the range is build as ub - lb.
        # For a constraint with an upper bound ('lower than'), the bound in rb
        # is used. For equality constraints, lb is used.
        
        rng      <- ct == "R"
        lbc      <- ct == "L"
        clb[lbc] <- ub[lbc]
        
        out <- cplexAPI::chgRhsCPLEX(lp@oobj@env, lp@oobj@lp,
                                     length(i), i-1, clb)

        if (sum(rng) > 0) {
            rngv <- abs(ub[rng] - lb[rng])
            out  <- cplexAPI::chgRngValCPLEX(lp@oobj@env, lp@oobj@lp,
                                   sum(rng), i[rng]-1, rngv)
        }
        
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setRhsZero", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        nrows  <- cplexAPI::getNumRowsCPLEX(lp@oobj@env, lp@oobj@lp)
        zeros  <- rep(0, nrows)
        indic  <- c(0:(nrows-1))
        outb   <- cplexAPI::chgRhsCPLEX(lp@oobj@env, lp@oobj@lp,
                                        nrows, indic, zeros)
        outt   <- cplexAPI::chgSenseCPLEX(lp@oobj@env, lp@oobj@lp,
                                          nrows, indic, rep("E", nrows))
        out <- c(outb, outt)
        # usable only for model creation!
        # ( Variable nrows has to be argument of setRhsZero()! )
        # out <- cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp, nrows,
        #                               rep(0, nrows), rep("E", nrows))

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsLowBnds", signature(lp = "optObj_cplexAPI", i = "numeric"),

    function(lp, i) {

        wrong_solver_msg(lp, "getRowsLowBnds", printOut = TRUE)
        return(FALSE)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsUppBnds", signature(lp = "optObj_cplexAPI", i = "numeric"),

    function(lp, i) {

        wrong_solver_msg(lp, "getRowsUppBnds", printOut = TRUE)
        return(FALSE)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeObjCoefs", signature(lp = "optObj_cplexAPI"),

    function(lp, j, obj_coef) {

        out <- cplexAPI::chgObjCPLEX(lp@oobj@env, lp@oobj@lp,
                                     length(j), j-1, obj_coef)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjCoefs", signature(lp = "optObj_cplexAPI", j = "numeric"),

    function(lp, j) {

        if (length(j) > 1) {
            b <- min(j) - 1
            e <- max(j) - 1
        }
        else {
            b <- j - 1
            e <- j - 1
        }
        out <- cplexAPI::getObjCPLEX(lp@oobj@env, lp@oobj@lp, b, e)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeMatrixRow", signature(lp = "optObj_cplexAPI"),

    function(lp, i, j, val) {

        cplexAPI::chgCoefListCPLEX(lp@oobj@env, lp@oobj@lp,
                                   length(val), rep(i-1, length(val)), j-1, val)

    }
)


#------------------------------------------------------------------------------#

#setMethod("loadLPprob", signature(lp = "optObj_cplexAPI"),
#
#    function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype,
#             lpdir = "max", rub = NULL, ctype = NULL,
#             cnames = NULL, rnames = NULL) {
#
#        stopifnot(is(mat, "Matrix"))
#
#        crtype <- sapply(rtype,
#                         function(x) switch(x,
#                                            "L" = { "G" },
#                                            "U" = { "L" },
#                                            "E" = { "E" },
#                                            "R" = { "R" },
#                                                  { "E" }))
#
#        # ranged constraints
#        if (is.null(rub)) {
#            crub <- NULL
#        }
#        else {
#            #rng        <- rtype == "R"
#            rng        <- rtype %in% "R"
#            crub       <- numeric(nRows)
#            crub[rng]  <- rlb[rng] - rlb[rng]
#            crub[!rng] <- 0
#        }
#
##
##        # problem type
##         ptype <- switch(lp@probType,
##             "lp"  = { CPXPROB_LP },
##             "mip" = { CPXPROB_MILP },
##                     { CPXPROB_LP }
##         )
##         cplexAPI::chgProbTypeCPLEX(lp@oobj@env, lp@oobj@lp, ptype)
#
#
##
##        # load problem
##        TMPmat <- as(mat, "CsparseMatrix")
##        cplexAPI::copyLpCPLEX(lp@oobj@env, lp@oobj@lp,
##                              nCols  = nCols,
##                              nRows  = nRows,
##                              lpdir  = ifelse(lpdir == "max",
##                                              cplexAPI::CPX_MAX,
##                                              cplexAPI::CPX_MIN),
##                              objf   = obj,
##                              rhs    = rlb,
##                              sense  = crtype,
##                              matbeg = TMPmat@p,
##                              matcnt = colSums(mat != 0),
##                              matind = TMPmat@i,
##                              matval = TMPmat@x,
##                              lb     = lb,
##                              ub     = ub,
##                              rngval = crub)
##
##        if (!is.null(ctype)) {
##            cplexAPI::chgColTypeCPLEX(lp@oobj@env, lp@oobj@lp,
##                                      ncols  = nCols,
##                                      ind    = c(1:nCols),
##                                      xctype = ctype)
##        }
##
#
#        # optimization direction
#        setObjDir(lp, lpdir = lpdir)
#
#        # constraints and right hand side
#        cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp,
#                               nrows  = nRows,
#                               rhs    = rlb,
#                               sense  = crtype,
#                               rngval = crub,
#                               rnames = rnames)
#
#        # variables, bounds and objective function
#        cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp,
#                               ncols  = nCols,
#                               obj    = obj,
#                               lb     = lb,
#                               ub     = ub,
#                               cnames = cnames)
#
#        # constraint matrix
#        TMPmat <- as(mat, "TsparseMatrix")
#        cplexAPI::chgCoefListCPLEX(lp@oobj@env, lp@oobj@lp,
#                                   nnz = length(TMPmat@x),
#                                   ia  = TMPmat@i,
#                                   ja  = TMPmat@j,
#                                   ra  = TMPmat@x)
#
#        if (!is.null(ctype)) {
#            cplexAPI::copyColTypeCPLEX(lp@oobj@env, lp@oobj@lp,
#                                       xctype = ctype)
#        }
#    }
#)


#------------------------------------------------------------------------------#

setMethod("loadLPprob", signature(lp = "optObj_cplexAPI"),

    function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype,
             lpdir = "max", rub = NULL, ctype = NULL,
             cnames = NULL, rnames = NULL, pname = NULL) {

        stopifnot(is(mat, "Matrix"))

        crtype <- sapply(rtype,
                         function(x) switch(EXPR = x,
                                            "L" = { "G" },
                                            "U" = { "L" },
                                            "D" = { "R" },
                                            "E" = { "E" },
                                                  { "E" }))

        # ranged constraints
        if (is.null(rub)) {
            crub <- NULL
            crlb <- rlb
        }
        else {
            # CPLEX only has a right-hand-side (rhs) and a so called range value
            # for reanged constraints. The value in rub is used to calculate the
            # range value (if required).
            # Range for constraint i is abs(rub[i] - rlb[i]) The interval for
            # constraint i then is [ rlb[i] , rlb[i] + range ] .
            # For constraints with an upper bound, the value in rub is copied
            # to rlb.
            
            stopifnot(length(rlb) == length(rub))
            rng        <- crtype == "R"
            crub       <- numeric(nRows)
            crub[rng]  <- abs(rub[rng] - rlb[rng])     # range
            
            ubc        <- crtype == "L"
            crlb       <- rlb
            crlb[ubc]  <- rub[ubc]
        }

#
#        # problem type
#         ptype <- switch(lp@probType,
#             "lp"  = { CPXPROB_LP },
#             "mip" = { CPXPROB_MILP },
#                     { CPXPROB_LP }
#         )
#         cplexAPI::chgProbTypeCPLEX(lp@oobj@env, lp@oobj@lp, ptype)


#
#        # load problem
#        TMPmat <- as(mat, "CsparseMatrix")
#        cplexAPI::copyLpCPLEX(lp@oobj@env, lp@oobj@lp,
#                              nCols  = nCols,
#                              nRows  = nRows,
#                              lpdir  = ifelse(lpdir == "max",
#                                              cplexAPI::CPX_MAX,
#                                              cplexAPI::CPX_MIN),
#                              objf   = obj,
#                              rhs    = rlb,
#                              sense  = crtype,
#                              matbeg = TMPmat@p,
#                              matcnt = colSums(mat != 0),
#                              matind = TMPmat@i,
#                              matval = TMPmat@x,
#                              lb     = lb,
#                              ub     = ub,
#                              rngval = crub)
#
#        if (!is.null(ctype)) {
#            cplexAPI::chgColTypeCPLEX(lp@oobj@env, lp@oobj@lp,
#                                      ncols  = nCols,
#                                      ind    = c(1:nCols),
#                                      xctype = ctype)
#        }
#

        # optimization direction
        setObjDir(lp, lpdir = lpdir)

        # constraints and right hand side
        cplexAPI::newRowsCPLEX(lp@oobj@env, lp@oobj@lp,
                               nrows  = nRows,
                               rhs    = crlb,
                               sense  = crtype,
                               rngval = crub,
                               rnames = rnames)

        # variables, bounds and objective function
        cplexAPI::newColsCPLEX(lp@oobj@env, lp@oobj@lp,
                               ncols  = nCols,
                               obj    = obj,
                               lb     = lb,
                               ub     = ub,
                               cnames = cnames)

        # constraint matrix
        TMPmat <- as(mat, "TsparseMatrix")
        cplexAPI::chgCoefListCPLEX(lp@oobj@env, lp@oobj@lp,
                                   nnz = length(TMPmat@x),
                                   ia  = TMPmat@i,
                                   ja  = TMPmat@j,
                                   ra  = TMPmat@x)

        if (!is.null(ctype)) {
            cplexAPI::copyColTypeCPLEX(lp@oobj@env, lp@oobj@lp,
                                       xctype = ctype)
        }

        # problem name
        if (!is.null(pname)) {
            cplexAPI::chgProbNameCPLEX(lp@oobj@env, lp@oobj@lp,
                                       probname = pname)
        }

    }
)


#------------------------------------------------------------------------------#

setMethod("loadQobj", signature(lp = "optObj_cplexAPI", mat = "Matrix"),

    function(lp, mat) {

        TMPmat <- as(mat, "CsparseMatrix")
        cplexAPI::copyQuadCPLEX(lp@oobj@env, lp@oobj@lp,
                                qmatbeg = TMPmat@p,
                                qmatcnt = colSums(mat != 0),
                                qmatind = TMPmat@i,
                                qmatval = TMPmat@x)

    }
)


#------------------------------------------------------------------------------#

setMethod("loadQobj", signature(lp = "optObj_cplexAPI", mat = "numeric"),

    function(lp, mat) {

        cplexAPI::copyQPsepCPLEX(lp@oobj@env, lp@oobj@lp, qsepvec = mat)

    }
)


#------------------------------------------------------------------------------#

setMethod("scaleProb", signature(lp = "optObj_cplexAPI"),

    function(lp, opt) {

        out <- cplexAPI::setIntParmCPLEX(lp@oobj@env,
                                         cplexAPI::CPX_PARAM_REDUCE,
                                         opt)
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("solveLp", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        out <- FALSE
        switch(EXPR = lp@method,
            "primopt" = {
                out <- cplexAPI::primoptCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            "dualopt" = {
                out <- cplexAPI::dualoptCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            "baropt" = {
                out <- cplexAPI::baroptCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            "hybbaropt" = {
                out <- cplexAPI::hybbaroptCPLEX(lp@oobj@env, lp@oobj@lp,
                                                method = 0)
            },
            "hybnetopt" = {
                out <- cplexAPI::hybnetoptCPLEX(lp@oobj@env, lp@oobj@lp,
                                      method = cplexAPI::CPX_ALG_PRIMAL)
            },
            "siftopt" = {
                out <- cplexAPI::siftoptCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            "mipopt" = {
                out <- cplexAPI::mipoptCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            "qpopt" = {
                out <- cplexAPI::qpoptCPLEX(lp@oobj@env, lp@oobj@lp)
            },
            {
                out <- cplexAPI::lpoptCPLEX(lp@oobj@env, lp@oobj@lp)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjVal", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        obj <- cplexAPI::getObjValCPLEX(lp@oobj@env, lp@oobj@lp)

        if (is(obj, "cplexError")) {
            if (probType(lp) == "mip") {
                out <- cplexAPI::getBestObjValCPLEX(lp@oobj@env, lp@oobj@lp)
            }
            else {
                out <- as.numeric(NA)
                #out <- 0
            }
        }
        else {
            out <- obj
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRedCosts", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        nr  <- cplexAPI::getNumRowsCPLEX(lp@oobj@env, lp@oobj@lp)
        out <- cplexAPI::getDjCPLEX(lp@oobj@env, lp@oobj@lp, 0, nr-1)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolStat", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        out <- cplexAPI::getStatCPLEX(lp@oobj@env, lp@oobj@lp)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getFluxDist", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        ncols <- cplexAPI::getNumColsCPLEX(lp@oobj@env, lp@oobj@lp)
        fluxd <- cplexAPI::getProbVarCPLEX(lp@oobj@env, lp@oobj@lp, 0, ncols-1)

        if (is(fluxd, "cplexError")) {
            out <- numeric(ncols)
        }
        else {
            out <- fluxd
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColPrim", signature(lp = "optObj_cplexAPI", j = "numeric"),

    function(lp, j) {

        out <- cplexAPI::getProbVarCPLEX(lp@oobj@env, lp@oobj@lp, j-1, j-1)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumNnz", signature(lp = "optObj_cplexAPI"),

    function(lp) {

        out <- cplexAPI::getNumNnzCPLEX(lp@oobj@env, lp@oobj@lp)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("writeProb", signature(lp = "optObj_cplexAPI", fname = "character"),

    function(lp, fname, ff = "lp") {

        tp  <- ifelse(is.null(ff), NULL, toupper(ff))
        fl  <- cplexAPI::writeProbCPLEX(lp@oobj@env, lp@oobj@lp,
                                        fname = fname, ftype = tp)
        out <- ifelse(fl == 0, TRUE, fl)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("readProb", signature(lp = "optObj_cplexAPI", fname = "character"),

    function(lp, fname, ff = "lp") {

        tp  <- ifelse(is.null(ff), NULL, toupper(ff))
        fl  <- cplexAPI::readCopyProbCPLEX(lp@oobj@env, lp@oobj@lp,
                                           fname = fname, ftype = tp)
        out <- ifelse(fl == 0, TRUE, fl)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("sensitivityAnalysis", signature(lp = "optObj_cplexAPI"),

    function(lp, ...) {

        # number of columns and rows
        nc <- cplexAPI::getNumColsCPLEX(lp@oobj@env, lp@oobj@lp)
        nr <- cplexAPI::getNumRowsCPLEX(lp@oobj@env, lp@oobj@lp)

        out <- vector(mode = "list", length = 3)
        names(out) <- c("bound", "obj", "rhs")
        
        out[["bound"]] <- cplexAPI::boundSaCPLEX(lp@oobj@env,
                                                 lp@oobj@lp, 0, nc-1)
        out[["obj"]]   <- cplexAPI::objSaCPLEX(lp@oobj@env,
                                               lp@oobj@lp, 0, nc-1)
        out[["rhs"]]   <- cplexAPI::rhsSaCPLEX(lp@oobj@env,
                                               lp@oobj@lp, 0, nr-1)

        return(out)
    }
)

#------------------------------------------------------------------------------#


setMethod("setRowsNames", signature(lp = "optObj_cplexAPI",
                                    i = "numeric", names = "character"),

    function(lp, i, names) {

        invisible(cplexAPI::chgRowNameCPLEX(lp@oobj@env, lp@oobj@lp,
                                            length(i), i-1, names))

    }
)


#------------------------------------------------------------------------------#

setMethod("setColsNames", signature(lp = "optObj_cplexAPI",
                                    j = "numeric", names = "character"),

    function(lp, j, names) {

        invisible(cplexAPI::chgColNameCPLEX(lp@oobj@env, lp@oobj@lp,
                                            length(j), j-1, names))

    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsNames", signature(lp = "optObj_cplexAPI", i = "numeric"),

    function(lp, i) {

        rn <- mapply(cplexAPI::getRowNameCPLEX, begin = i-1, end = i-1,
                     MoreArgs = list(env = lp@oobj@env, lp = lp@oobj@lp),
                     SIMPLIFY = TRUE, USE.NAMES = FALSE)
        return(unlist(rn))

    }
)


#------------------------------------------------------------------------------#

setMethod("getColsNames", signature(lp = "optObj_cplexAPI", j = "numeric"),

    function(lp, j) {

        cn <- mapply(cplexAPI::getColNameCPLEX, begin = j-1, end = j-1,
                     MoreArgs = list(env = lp@oobj@env, lp = lp@oobj@lp),
                     SIMPLIFY = TRUE, USE.NAMES = FALSE)
        return(unlist(cn))

    }
)


#------------------------------------------------------------------------------#

