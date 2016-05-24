#  optObj_glpkAPIClass.R
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
#                  definition of the class optObj_glpkAPI                      #
#------------------------------------------------------------------------------#

setClass(Class = "optObj_glpkAPI", contains = "optObj")


#------------------------------------------------------------------------------#
#                                  methods                                     #
#------------------------------------------------------------------------------#

setMethod("delProb", signature(lp = "optObj_glpkAPI"),

    function(lp, ...) {

        glpkAPI::delProbGLPK(lp@oobj)

    }
)


#------------------------------------------------------------------------------#

setMethod("initProb", signature(lp = "optObj_glpkAPI"),

    function(lp, to = FALSE, ...) {

        lp@oobj <- glpkAPI::initProbGLPK()

        if (is.null(to)) {
            too <- FALSE
        }
        else {
            too <- to
        }
        
        if (isTRUE(too)) {
            glpkAPI::termOutGLPK(glpkAPI::GLP_ON)
        }
        else {
            glpkAPI::termOutGLPK(glpkAPI::GLP_OFF)
        }

        return(lp)
    }
)


#------------------------------------------------------------------------------#

setMethod("backupProb", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        out <- FALSE
        np  <- FALSE

        # reset parameters!!! Parameters are reset to default when
        # doing initProbGLPK() !!!
        np <- glpkAPI::initProbGLPK()
        glpkAPI::copyProbGLPK(lp@oobj, np)

        # create new optObj object
        if (!identical(np, FALSE)) {
            out <- new("optObj_glpkAPI", lp@solver, lp@method, lp@probType)
            out@oobj <- np
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setSolverParm", signature(lp = "optObj_glpkAPI"),

    function(lp, solverParm) {

        if ( ! ((is.data.frame(solverParm)) || (is.list(solverParm))) ) {
            stop(sQuote(solverParm), " must be list or data.frame")
        }

        if (any(is.na(solverParm))) {
            stop(sQuote(solverParm), " contains NA values")
        }

        parm <- sapply(names(solverParm), function(x) eval(parse(text = x)))
        val  <- unlist(solverParm)

        switch(EXPR=lp@method,
            "interior" = {
                glpkAPI::setInteriorParmGLPK(parm, val)
            },
            "mip" = {
                glpkAPI::setMIPParmGLPK(parm, val)
            },
            {
                glpkAPI::setSimplexParmGLPK(parm, val)
            }
        )
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolverParm", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        out <- FALSE

        switch(EXPR=lp@method,
            "interior" = {
                out <- glpkAPI::getInteriorParmGLPK()
            },
            "mip" = {
                out <- glpkAPI::getMIPParmGLPK()
            },
            {
                out <- glpkAPI::getSimplexParmGLPK()
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_glpkAPI", lpdir = "character"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == "max", glpkAPI::GLP_MAX, glpkAPI::GLP_MIN)
        glpkAPI::setObjDirGLPK(lp@oobj, dr)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_glpkAPI", lpdir = "integer"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == glpkAPI::GLP_MAX,
                     glpkAPI::GLP_MAX,
                     glpkAPI::GLP_MIN)
        glpkAPI::setObjDirGLPK(lp@oobj, dr)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_glpkAPI", lpdir = "numeric"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == -1,
                     glpkAPI::GLP_MAX,
                     glpkAPI::GLP_MIN)
        glpkAPI::setObjDirGLPK(lp@oobj, dr)
    }
)

#------------------------------------------------------------------------------#

setMethod("getObjDir", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        dr <- glpkAPI::getObjDirGLPK(lp@oobj)
        if (dr == glpkAPI::GLP_MAX) {
            out <- "max"
        }
        else if (dr == glpkAPI::GLP_MIN) {
            out <- "min"
        }
        else {
            out <- FALSE
        }
        
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRows", signature(lp = "optObj_glpkAPI", nrows = "numeric"),

    function(lp, nrows) {

        out <- glpkAPI::addRowsGLPK(lp@oobj, nrows)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addCols", signature(lp = "optObj_glpkAPI", ncols = "numeric"),

    function(lp, ncols) {

        out <- glpkAPI::addColsGLPK(lp@oobj, ncols)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsCols", signature(lp = "optObj_glpkAPI",
                                   nrows = "numeric", ncols = "numeric"),

    function(lp, nrows, ncols) {

        outi <- glpkAPI::addRowsGLPK(lp@oobj, nrows)
        outj <- glpkAPI::addColsGLPK(lp@oobj, ncols)
        out  <- all(outi, outj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumRows", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        out <- glpkAPI::getNumRowsGLPK(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumCols", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        out <- glpkAPI::getNumColsGLPK(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addColsToProb", signature(lp = "optObj_glpkAPI"),

    # j: vector containing the new column indices (must be ascending)
    # rind: list, containing the row indices of the new nz elements
    # nzval: list, containing the new nz elements
    #
    # j, obj, lb, rind and nzval must have the same length

    function(lp, j, obj, lb, ub, rind, nzval) {

        ord <- glpkAPI::addColsGLPK(lp@oobj, length(j))
        for (k in seq(along = j)) {
            glpkAPI::setMatColGLPK(lp@oobj, j[k],
                                   length(rind[[k]]),
                                   rind[[k]], nzval[[k]])
        }
        out <- glpkAPI::setColsBndsGLPK(lp@oobj, j, lb, ub)
        out <- glpkAPI::setObjCoefsGLPK(lp@oobj, j, obj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsToProb", signature(lp = "optObj_glpkAPI"),

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

        ord <- glpkAPI::addRowsGLPK(lp@oobj, length(i))
        gtype = integer(length(type))
        for (k in seq(along = i)) {
            gtype[k] <- switch(EXPR = type[k],
                               "F" = { glpkAPI::GLP_FR },
                               "L" = { glpkAPI::GLP_LO },
                               "U" = { glpkAPI::GLP_UP },
                               "D" = { glpkAPI::GLP_DB },
                               "E" = { glpkAPI::GLP_FX },
                                     { glpkAPI::GLP_FX }
            )
            glpkAPI::setMatRowGLPK(lp@oobj, i[k],
                                   length(cind[[k]]),
                                   cind[[k]], nzval[[k]])
        }
        stopifnot(length(lb) == length(ub))
        out <- glpkAPI::setRowsBndsGLPK(lp@oobj, i, lb, ub, gtype)

        # row names
        if (!is.null(rnames)) {
            glpkAPI::setRowsNamesGLPK(lp@oobj, i = i, rnames = rnames)
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBnds", signature(lp = "optObj_glpkAPI"),

    function(lp, j, lb, ub) {

        glpkAPI::setColsBndsGLPK(lp@oobj, j, lb, ub)

    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBndsObjCoefs", signature(lp = "optObj_glpkAPI"),

    function(lp, j, lb, ub, obj_coef) {

        glpkAPI::setColsBndsObjCoefsGLPK(lp@oobj, j, lb, ub, obj_coef)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsLowBnds", signature(lp = "optObj_glpkAPI", j = "numeric"),

    function(lp, j) {

        out <- glpkAPI::getColsLowBndsGLPK(lp@oobj, j)
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsUppBnds", signature(lp = "optObj_glpkAPI", j = "numeric"),

    function(lp, j) {

        out <- glpkAPI::getColsUppBndsGLPK(lp@oobj, j)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeRowsBnds", signature(lp = "optObj_glpkAPI"),

    function(lp, i, lb, ub) {

        ct <- glpkAPI::getRowsTypesGLPK(lp@oobj, i)
        glpkAPI::setRowsBndsGLPK(lp@oobj, i, lb, ub, type = ct)

    }
)


#------------------------------------------------------------------------------#

setMethod("setRhsZero", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        glpkAPI::setRhsZeroGLPK(lp@oobj)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsLowBnds", signature(lp = "optObj_glpkAPI", i = "numeric"),

    function(lp, i) {

        out <- glpkAPI::getRowsLowBndsGLPK(lp@oobj, i)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsUppBnds", signature(lp = "optObj_glpkAPI", i = "numeric"),

    function(lp, i) {

        out <- glpkAPI::getRowsUppBndsGLPK(lp@oobj, i)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeObjCoefs", signature(lp = "optObj_glpkAPI"),

    function(lp, j, obj_coef) {

        glpkAPI::setObjCoefsGLPK(lp@oobj, j, obj_coef)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjCoefs", signature(lp = "optObj_glpkAPI", j = "numeric"),

    function(lp, j) {

        out <- glpkAPI::getObjCoefsGLPK(lp@oobj, j)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeMatrixRow", signature(lp = "optObj_glpkAPI"),

    function(lp, i, j, val) {

        glpkAPI::setMatRowGLPK(lp@oobj, i, length(val), j, val)

    }
)


#------------------------------------------------------------------------------#

#setMethod("loadLPprob", signature(lp = "optObj_glpkAPI"),
#
#    function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype,
#             lpdir = "max", rub = NULL, ctype = NULL,
#             cnames = NULL, rnames = NULL) {
#
#        stopifnot(is(mat, "Matrix"))
#
#        crtype <- sapply(rtype,
#                         function(x) switch(x,
#                                            "F" = { glpkAPI::GLP_FR },
#                                            "L" = { glpkAPI::GLP_LO },
#                                            "U" = { glpkAPI::GLP_UP },
#                                            "D" = { glpkAPI::GLP_DB },
#                                            "E" = { glpkAPI::GLP_FX },
#                                                  { glpkAPI::GLP_FX }))
#
#        # optimization direction
#        setObjDir(lp, lpdir = lpdir)
#        
#        # problem dimensions
#        glpkAPI::addColsGLPK(lp@oobj, ncols = nCols)
#        glpkAPI::addRowsGLPK(lp@oobj, nrows = nRows)
#
#        # constraint matrix
#        TMPmat <- as(mat, "TsparseMatrix")
#        glpkAPI::loadMatrixGLPK(lp@oobj,
#                                ne = length(TMPmat@x),
#                                ia = TMPmat@i + 1,
#                                ja = TMPmat@j + 1,
#                                ra = TMPmat@x)
#
#        # column (variable) bounds and objective function
#        glpkAPI::setColsBndsObjCoefsGLPK(lp@oobj,
#                                         j = c(1:nCols),
#                                         lb = lb,
#                                         ub = ub,
#                                         obj_coef = obj)
#
#        # variable type
#        if (!is.null(ctype)) {
#            cctype <- sapply(ctype,
#                             function(x) switch(x,
#                                                "C" = { glpkAPI::GLP_CV },
#                                                "I" = { glpkAPI::GLP_IV },
#                                                "B" = { glpkAPI::GLP_BV },
#                                                      { glpkAPI::GLP_CV }))
#
#            glpkAPI::setColsKindGLPK(lp@oobj, j = c(1:nCols), kind = cctype)
#        }
#
#        # right hand side
#        glpkAPI::setRowsBndsGLPK(lp@oobj,
#                                 i = c(1:nRows),
#                                 lb = rlb,
#                                 ub = rub,
#                                 type = crtype)
#
#        # row names
#        if (!is.null(rnames)) {
#            glpkAPI::setRowsNamesGLPK(lp@oobj, i = c(1:nRows), rnames = rnames)
#        }
#
#        # column names
#        if (!is.null(cnames)) {
#            glpkAPI::setColsNamesGLPK(lp@oobj, j = c(1:nCols), cnames = cnames)
#        }
#
#        if (!is.null(rnames) || !is.null(cnames)) {
#            glpkAPI::createIndexGLPK(lp@oobj)
#        }
#    }
#)


#------------------------------------------------------------------------------#

setMethod("loadLPprob", signature(lp = "optObj_glpkAPI"),

    function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype,
             lpdir = "max", rub = NULL, ctype = NULL,
             cnames = NULL, rnames = NULL, pname = NULL) {

        stopifnot(is(mat, "Matrix"))

        crtype <- sapply(rtype,
                         function(x) switch(EXPR = x,
                                            "F" = { glpkAPI::GLP_FR },
                                            "L" = { glpkAPI::GLP_LO },
                                            "U" = { glpkAPI::GLP_UP },
                                            "D" = { glpkAPI::GLP_DB },
                                            "E" = { glpkAPI::GLP_FX },
                                                  { glpkAPI::GLP_FX }))

        
        # optimization direction
        setObjDir(lp, lpdir = lpdir)
        
        # problem dimensions
        glpkAPI::addColsGLPK(lp@oobj, ncols = nCols)
        glpkAPI::addRowsGLPK(lp@oobj, nrows = nRows)

        # constraint matrix
        TMPmat <- as(mat, "TsparseMatrix")
        glpkAPI::loadMatrixGLPK(lp@oobj,
                                ne = length(TMPmat@x),
                                ia = TMPmat@i + 1,
                                ja = TMPmat@j + 1,
                                ra = TMPmat@x)

        # column (variable) bounds and objective function
        glpkAPI::setColsBndsObjCoefsGLPK(lp@oobj,
                                         j = c(1:nCols),
                                         lb = lb,
                                         ub = ub,
                                         obj_coef = obj)

        # variable type
        if (!is.null(ctype)) {
            cctype <- sapply(ctype,
                             function(x) switch(EXPR = x,
                                                "C" = { glpkAPI::GLP_CV },
                                                "I" = { glpkAPI::GLP_IV },
                                                "B" = { glpkAPI::GLP_BV },
                                                      { glpkAPI::GLP_CV }))

            glpkAPI::setColsKindGLPK(lp@oobj, j = c(1:nCols), kind = cctype)
        }

        # right hand side
        if (is.null(rub)) {
            # The values in rlb will be copied to rub. GLPK ignores rlb and rub,
            # depending on the constraint type (e.g. an upper bound, if the
            # constraint type says, it has a lower bound):
            # Constraint type "L": ignore rub
            # Constraint type "U": ignore rlb
            # Constraint type "E": ignore rub
            # Constraint type "F": ignore rlb and rub

            crub <- rlb
        }
        else {
            crub <- rub
        }
        stopifnot(length(rlb) == length(crub))
        glpkAPI::setRowsBndsGLPK(lp@oobj,
                                 i = c(1:nRows),
                                 lb = rlb,
                                 ub = crub,
                                 type = crtype)

        # row names
        if (!is.null(rnames)) {
            glpkAPI::setRowsNamesGLPK(lp@oobj, i = c(1:nRows), rnames = rnames)
        }

        # column names
        if (!is.null(cnames)) {
            glpkAPI::setColsNamesGLPK(lp@oobj, j = c(1:nCols), cnames = cnames)
        }

        # problem name
        if (!is.null(pname)) {
            glpkAPI::setProbNameGLPK(lp@oobj, pname = pname)
        }

        if (!is.null(rnames) || !is.null(cnames) || !is.null(pname)) {
            glpkAPI::createIndexGLPK(lp@oobj)
        }
    }
)


#------------------------------------------------------------------------------#

setMethod("scaleProb", signature(lp = "optObj_glpkAPI"),

    function(lp, opt) {

        # check if tryCatch works here!!
        glpkAPI::scaleProbGLPK(lp@oobj, opt)

    }
)


#------------------------------------------------------------------------------#

setMethod("solveLp", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        out          <- FALSE
#        if (glpkAPI::bfExistsGLPK(lp@oobj) != 0) {
#            if (glpkAPI::bfUpdatedGLPK(lp@oobj) != 0) {
#                basis <- glpkAPI::factorizeGLPK(lp@oobj)
#                #print(basis)
#            }
#        }
        switch(EXPR = lp@method,
            "interior" = {
                out <- glpkAPI::solveInteriorGLPK(lp@oobj)
            },
            "exact" = {
                out <- glpkAPI::solveSimplexExactGLPK(lp@oobj)
            },
            "mip" = {
                out <- glpkAPI::solveMIPGLPK(lp@oobj)
            },
            {
                out <- glpkAPI::solveSimplexGLPK(lp@oobj)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjVal", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        switch(EXPR = lp@method,
            "interior" = {
                out <- glpkAPI::getObjValIptGLPK(lp@oobj)
            },
            "mip" = {
                out <- glpkAPI::mipObjValGLPK(lp@oobj)
            },
            {
                out <- glpkAPI::getObjValGLPK(lp@oobj)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRedCosts", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        if (lp@method == "interior") {
            out <- glpkAPI::getColsDualIptGLPK(lp@oobj)
        }
        else {
            out <- glpkAPI::getColsDualGLPK(lp@oobj)
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolStat", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        switch(EXPR = lp@method,
            "interior" = {
                out <- glpkAPI::getSolStatIptGLPK(lp@oobj)
            },
            "mip" = {
                out <- glpkAPI::mipStatusGLPK(lp@oobj)
            },
            {
                out <- glpkAPI::getSolStatGLPK(lp@oobj)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getFluxDist", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        switch(EXPR = lp@method,
            "interior" = {
                out <- glpkAPI::getColsPrimIptGLPK(lp@oobj)
            },
            "mip" = {
                out <- glpkAPI::mipColsValGLPK(lp@oobj)
            },
            {
                out <- glpkAPI::getColsPrimGLPK(lp@oobj)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColPrim", signature(lp = "optObj_glpkAPI", j = "numeric"),

    function(lp, j) {

        switch(EXPR = lp@method,
            "interior" = {
                out <- glpkAPI::getColPrimIptGLPK(lp@oobj, j)
            },
            "mip" = {
                out <- glpkAPI::mipColValGLPK(lp@oobj, j)
            },
            {
                out <- glpkAPI::getColPrimGLPK(lp@oobj, j)
            }
        )

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumNnz", signature(lp = "optObj_glpkAPI"),

    function(lp) {

        out <- glpkAPI::getNumNnzGLPK(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("writeProb", signature(lp = "optObj_glpkAPI", fname = "character"),

    function(lp, fname, ff = "lp", ...) {

        switch(EXPR = ff,
            "lp"  = {
                fl <- glpkAPI::writeLPGLPK(lp@oobj, fname = fname)
            },
            "mps" = {
                fl <- glpkAPI::writeMPSGLPK(lp@oobj, fname = fname, ...)
            },
            "glpk" = {
                fl <- glpkAPI::writeProbGLPK(lp@oobj, fname = fname)
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

setMethod("readProb", signature(lp = "optObj_glpkAPI", fname = "character"),

    function(lp, fname, ff = "lp", ...) {

        switch(EXPR = ff,
            "lp"  = {
                fl <- glpkAPI::readLPGLPK(lp@oobj, fname = fname)
            },
            "mps" = {
                fl <- glpkAPI::readMPSGLPK(lp@oobj, fname = fname, ...)
            },
            "glpk" = {
                fl <- glpkAPI::readProbGLPK(lp@oobj, fname = fname)
            },
            {
                message("wrong format!")
                fl <- 1
            }
        )
        out <- ifelse(fl == 0, TRUE, fl)

        return(fl)
    }
)


#------------------------------------------------------------------------------#

setMethod("sensitivityAnalysis", signature(lp = "optObj_glpkAPI"),

    function(lp, ...) {

        out <- glpkAPI::printRangesGLPK(lp@oobj, ...)
        if (out == 0) {
            message("wrote the file 'sar.txt'")
        }
        else {
            warning("sensitivity analysis failed")
        }

        return(out)
    }
)

#------------------------------------------------------------------------------#


setMethod("setRowsNames", signature(lp = "optObj_glpkAPI",
                                    i = "numeric", names = "character"),

    function(lp, i, names) {

        invisible(glpkAPI::setRowsNamesGLPK(lp@oobj, i, names))

    }
)


#------------------------------------------------------------------------------#

setMethod("setColsNames", signature(lp = "optObj_glpkAPI",
                                    j = "numeric", names = "character"),

    function(lp, j, names) {

        invisible(glpkAPI::setColsNamesGLPK(lp@oobj, j, names))

    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsNames", signature(lp = "optObj_glpkAPI", i = "numeric"),

    function(lp, i) {

        rn <- mapply(glpkAPI::getRowNameGLPK, i, MoreArgs = list(lp = lp@oobj),
                     SIMPLIFY = TRUE, USE.NAMES = FALSE)
        return(unlist(rn))

    }
)


#------------------------------------------------------------------------------#

setMethod("getColsNames", signature(lp = "optObj_glpkAPI", j = "numeric"),

    function(lp, j) {

        cn <- mapply(glpkAPI::getColNameGLPK, j, MoreArgs = list(lp = lp@oobj),
                     SIMPLIFY = TRUE, USE.NAMES = FALSE)
        return(unlist(cn))

    }
)


#------------------------------------------------------------------------------#

