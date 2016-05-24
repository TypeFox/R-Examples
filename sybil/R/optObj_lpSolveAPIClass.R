#  optObj_lpSolveAPIClass.R
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
#                 definition of the class optObj_lpSolveAPI                    #
#------------------------------------------------------------------------------#

setClass(Class = "optObj_lpSolveAPI", contains = "optObj")


#------------------------------------------------------------------------------#
#                                  methods                                     #
#------------------------------------------------------------------------------#

setMethod("delProb", signature(lp = "optObj_lpSolveAPI"),

    function(lp, ...) {

        #finalizeLpSolveProb(lp@oobj)
        #lpSolveAPI::delete.lp(lp@oobj)
        #lp@oobj <- NULL
        lp <- new("optObj_lpSolveAPI")

        return(lp)

    }
)


#------------------------------------------------------------------------------#

setMethod("initProb", signature(lp = "optObj_lpSolveAPI"),

    function(lp, to = NULL, nrows = 0, ncols = 0) {

        lp@oobj <- lpSolveAPI::make.lp(nrow = nrows, ncol = ncols)
        #reg.finalizer(lp@oobj, finalizeLpSolveProb, TRUE)

        if (is.null(to)) {
            setSolverParm(lp, list(verbose = "neutral"))
        }
        else {
            setSolverParm(lp, list(verbose = as.character(to)))
        }

        return(lp)
    }
)


#------------------------------------------------------------------------------#

setMethod("backupProb", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- FALSE
        np  <- FALSE

        fname <- tempfile(pattern = "LPSOLVE_PROB_", fileext = ".tmp")
        lpSolveAPI::write.lp(lp@oobj, filename = fname, type = "lp")
        if (isTRUE(file.exists(fname))) {
            np <- lpSolveAPI::read.lp(fname, type = "lp")
            unlink(fname)
        }
        else {
            stop("cannot read model")
        }

        # create new optObj object
        if (!identical(np, FALSE)) {
            out <- new("optObj_lpSolveAPI", lp@solver, lp@method, lp@probType)
            out@oobj <- np
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setSolverParm", signature(lp = "optObj_lpSolveAPI"),

    function(lp, solverParm) {

        if ( ! ((is.data.frame(solverParm)) || (is.list(solverParm))) ) {
            stop(sQuote(solverParm), " must be list or data.frame")
        }

        if (any(is.na(solverParm))) {
            stop(sQuote(solverParm), " contains NA values")
        }

        pname <- names(solverParm)
        for (i in seq(along = solverParm)) {
            command <- paste("lpSolveAPI::lp.control(lp@oobj, ",
                             pname[i], "='" , solverParm[[i]], "')", sep = "")
            #print(command)
            eval(parse(text = command))
        }
        #print(lp.control(lp@oobj))

    }
)


#------------------------------------------------------------------------------#

setMethod("getSolverParm", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        lpSolveAPI::lp.control(lp@oobj)

    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_lpSolveAPI", lpdir = "character"),

    function(lp, lpdir) {

        lpSolveAPI::lp.control(lp@oobj, sense = lpdir)

    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_lpSolveAPI", lpdir = "numeric"),

    function(lp, lpdir) {

        dr <- ifelse(lpdir == -1, "max", "min")
        lpSolveAPI::lp.control(lp@oobj, sense = dr)

    }
)

#------------------------------------------------------------------------------#

setMethod("getObjDir", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- lpSolveAPI::lp.control(lp@oobj)[["sense"]]

        return(out)

    }
)


#------------------------------------------------------------------------------#

setMethod("addRows", signature(lp = "optObj_lpSolveAPI", nrows = "numeric"),

    function(lp, nrows) {

        ncols <- dim(lp@oobj)[2]
        out <- lpSolveAPI::resize.lp(lp@oobj, nrows, ncols)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addCols", signature(lp = "optObj_lpSolveAPI", ncols = "numeric"),

    function(lp, ncols) {

        nrows <- dim(lp@oobj)[1]
        out <- lpSolveAPI::resize.lp(lp@oobj, nrows, ncols)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsCols", signature(lp = "optObj_lpSolveAPI",
                                   nrows = "numeric", ncols = "numeric"),

    function(lp, nrows, ncols) {

        out <- lpSolveAPI::resize.lp(lp@oobj, nrows, ncols)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumRows", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- dim(lp@oobj)[1]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumCols", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- dim(lp@oobj)[2]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addColsToProb", signature(lp = "optObj_lpSolveAPI"),

    # j: vector containing the new column indices (must be ascending)
    # rind: list, containing the row indices of the new nz elements
    # nzval: list, containing the new nz elements
    #
    # j, obj, lb, rind and nzval must have the same length

    function(lp, j, obj, lb, ub, rind, nzval) {

        nrc <- dim(lp@oobj)
        obc <- getObjCoefs(lp, 1:nrc[2])

        ord <- lpSolveAPI::resize.lp(lp@oobj, nrc[1], (nrc[2]+length(j)))
        for (k in seq(along = j)) {
            lpSolveAPI::add.column(lp@oobj, nzval[[k]], rind[[k]])
        }
        out <- lpSolveAPI::set.bounds(lp@oobj, lb, ub, j)
        out <- lpSolveAPI::set.objfn(lp@oobj, c(obc, obj))

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsToProb", signature(lp = "optObj_lpSolveAPI"),

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
        for (k in seq(along = i)) {
            ltype <- switch(EXPR = type[k],
                            "L" = { 2 },
                            "U" = { 1 },
                            "D" = { 1 },
                            "E" = { 3 },
                                  { 3 }
            )


            if (type[k] == "D") {
                out <- lpSolveAPI::add.constraint(lp@oobj, nzval[[k]],
                                                 ltype, ub[k], cind[[k]], lb[k])
            }
            else {
                clb <- ifelse(type[k] == "U", ub[k], lb[k])
                out <- lpSolveAPI::add.constraint(lp@oobj, nzval[[k]],
                                                  ltype, clb, cind[[k]])
            }
        }

        # row names
        if (!is.null(rnames)) {
            rrnames <- sub("(", "_", rnames, fixed = TRUE)
            rrnames <- sub(")", "_", rrnames, fixed = TRUE)
            rownames(lp@oobj)[i] <- rrnames
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBnds", signature(lp = "optObj_lpSolveAPI"),

    function(lp, j, lb, ub) {

        out <- lpSolveAPI::set.bounds(lp@oobj, lb, ub, j)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBndsObjCoefs", signature(lp = "optObj_lpSolveAPI"),

    function(lp, j, lb, ub, obj_coef) {

        outb <- lpSolveAPI::set.bounds(lp@oobj, lb, ub, j)
        outo <- lpSolveAPI::set.objfn(lp@oobj, obj_coef, j)
        out  <- c(outb, outo)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsLowBnds", signature(lp = "optObj_lpSolveAPI", j = "numeric"),

    function(lp, j) {

        out <- lpSolveAPI::get.bounds(lp@oobj, j)[["lower"]]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsUppBnds", signature(lp = "optObj_lpSolveAPI", j = "numeric"),

    function(lp, j) {

        out <- lpSolveAPI::get.bounds(lp@oobj, j)[["upper"]]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeRowsBnds", signature(lp = "optObj_lpSolveAPI"),

    function(lp, i, lb, ub) {

        out <- lpSolveAPI::set.constr.value(lp@oobj, rhs = lb, lhs = ub,
                                                                constraints = i)
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setRhsZero", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        nrows <- dim(lp@oobj)[1]
        outb <- lpSolveAPI::set.constr.value(lp@oobj,
                                             rhs = rep(0, nrows),
                                             lhs = NULL,
                                             constraints = c(1:nrows))
        outt <- lpSolveAPI::set.constr.type(lp@oobj,
                                            rep(3, nrows), c(1:nrows))
        out <- c(outb, outt)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsLowBnds", signature(lp = "optObj_lpSolveAPI", i = "numeric"),

    function(lp, i) {

        out <- lpSolveAPI::get.constr.value(lp@oobj,
                                            side = "lhs",
                                            constraints = i)
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsUppBnds", signature(lp = "optObj_lpSolveAPI", i = "numeric"),

    function(lp, i) {

        out <- lpSolveAPI::get.constr.value(lp@oobj,
                                            side = "rhs",
                                            constraints = i)
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeObjCoefs", signature(lp = "optObj_lpSolveAPI"),

    function(lp, j, obj_coef) {

        out <- lpSolveAPI::set.objfn(lp@oobj, obj_coef, j)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjCoefs", signature(lp = "optObj_lpSolveAPI", j = "numeric"),

    function(lp, j) {

        #wrong_solver_msg(lp, "getObjCoefs", printOut = TRUE)
        out <- numeric(length(j))
        for (i in seq(along = j)) {
            out[i] <- lpSolveAPI::get.column(lp@oobj, j[i])$column[1]
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeMatrixRow", signature(lp = "optObj_lpSolveAPI"),

    function(lp, i, j, val) {

        lpSolveAPI::set.row(lp@oobj, row = i, xt = val, indices = j)

    }
)


#------------------------------------------------------------------------------#

setMethod("loadLPprob", signature(lp = "optObj_lpSolveAPI"),

    function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype,
             lpdir = "max", rub = NULL, ctype = NULL,
             cnames = NULL, rnames = NULL, pname = NULL) {

        stopifnot(is(mat, "Matrix"))

        # problem size is already specified in optObj_lpSolveAPI constructor

        crtype <- sapply(rtype,
                         function(x) switch(EXPR = x,
                                            "L" = { 2 },
                                            "U" = { 1 },
                                            "D" = { 1 },
                                            "E" = { 3 },
                                                  { 3 }))

        # constraints with double bound
        if (is.null(rub)) {
            crub <- NULL
            crlb <- rlb
        }
        else {
            # Package lpSolveAPI uses rhs and lhs for constraints, but -- as far
            # as I can see -- lhs does not work.
            # If there is a maximization problem, the values from rub for double
            # bounded constraints and constraints with upper bounds are
            # exchanged for the corresponding values in rlb.
            # If there is a minimization problem, the values from rub
            # constraints with only upper bounds are exchanged for the
            # corresponding values in rlb. For constraints with double bounds,
            # the constraint type is changed to a constraint with lower bound
            # (the upper bound should not be reached, because of the
            #  minimization).
            # Now, rlb is used as rhs and rub as rls for lpSolveAPI.
            
            stopifnot(length(rlb) == length(rub))
            crub       <- rub
            crlb       <- rlb

            # for a max problem: the constraining bound is the upper bound
            # for a min problem: the constraining bound is the lower bound

            if (lpdir == "max") {
                dbc    <- crtype == 1      # double bound and upper bound
            }
            else {
                dbc    <- rtype == "U"     # only upper bounds
                dbl    <- rtype == "D"
                crtype[dbl] <- 2
            }
            tmp        <- crlb[dbc]
            crlb[dbc]  <- crub[dbc]
            crub[dbc]  <- tmp
        }

        # optimization direction
        lpSolveAPI::lp.control(lp@oobj, sense = lpdir)

        # constraint matrix
        loadMatrixPerColumnLPSOLVE(lp@oobj, constMat = as(mat, "CsparseMatrix"))

        # column (variable) bounds
        lpSolveAPI::set.bounds(lp@oobj,
                               lower   = lb,
                               upper   = ub,
                               columns = c(1:nCols))

        # objective function
        lpSolveAPI::set.objfn(lp@oobj, obj = obj, indices = c(1:nCols))

        # constraints (right hand side)
        lpSolveAPI::set.constr.value(lp@oobj,
                                     rhs = crlb,
                                     lhs = crub,
                                     constraints = c(1:nRows))

        # type of constraint
        lpSolveAPI::set.constr.type(lp@oobj,
                                    type = crtype,
                                    constraints = c(1:nRows))

        # variable type
        if (!is.null(ctype)) {
            cont <- which(ctype %in% "C")
            int  <- which(ctype %in% "I")
            bin  <- which(ctype %in% "B")
            
            if (length(cont) > 0) {
                lpSolveAPI::set.type(lp@oobj, columns = cont, type = "real")
            }

            if (length(int) > 0) {
                lpSolveAPI::set.type(lp@oobj, columns = int, type = "integer")
            }

            if (length(bin) > 0) {
                lpSolveAPI::set.type(lp@oobj, columns = bin, type = "binary")
            }
        }

        # row and column names
        if ( (!is.null(rnames)) && (!is.null(cnames)) ) {
            rrnames <- sub("(", "_", rnames, fixed = TRUE)
            rrnames <- sub(")", "_", rrnames, fixed = TRUE)
            ccnames <- sub("(", "_", cnames, fixed = TRUE)
            ccnames <- sub(")", "_", ccnames, fixed = TRUE)
            dimnames(lp@oobj) <- list(rrnames, ccnames)
        }

        # problem name
        if (!is.null(pname)) {
            ppname <- sub("(", "_", pname, fixed = TRUE)
            ppname <- sub(")", "_", ppname, fixed = TRUE)
            lpSolveAPI::name.lp(lp@oobj, ppname)
        }
    }
)


#------------------------------------------------------------------------------#

setMethod("scaleProb", signature(lp = "optObj_lpSolveAPI"),

    function(lp, opt) {

        invisible(lpSolveAPI::lp.control(lp@oobj, scaling = opt))

    }
)


#------------------------------------------------------------------------------#

setMethod("solveLp", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- solve(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjVal", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- lpSolveAPI::get.objective(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRedCosts", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- lpSolveAPI::get.dual.solution(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolStat", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- NA
        
        wrong_solver_msg(lp, "getSolStat", printOut = FALSE)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getFluxDist", signature(lp = "optObj_lpSolveAPI"),

    function(lp) {

        out <- lpSolveAPI::get.variables(lp@oobj)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColPrim", signature(lp = "optObj_lpSolveAPI", j = "numeric"),

    function(lp, j) {

        out <- lpSolveAPI::get.variables(lp@oobj)[j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("writeProb", signature(lp = "optObj_lpSolveAPI", fname = "character"),

    function(lp, fname, ff = "lp", ...) {

        out <- lpSolveAPI::write.lp(lp@oobj, filename = fname, type = ff, ...)

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("readProb", signature(lp = "optObj_lpSolveAPI", fname = "character"),

    function(lp, fname, ff = "lp", ...) {

        lp@oobj <- lpSolveAPI::read.lp(filename = fname, type = ff, ...)

        return(lp)
    }
)


#------------------------------------------------------------------------------#

setMethod("setRowsNames", signature(lp = "optObj_lpSolveAPI",
                                    i = "numeric", names = "character"),

    function(lp, i, names) {

        invisible(dimnames(lp@oobj)[[1]][i] <- names)

    }
)


#------------------------------------------------------------------------------#

setMethod("setColsNames", signature(lp = "optObj_lpSolveAPI",
                                    j = "numeric", names = "character"),

    function(lp, j, names) {

        invisible(dimnames(lp@oobj)[[2]][j] <- names)

    }
)


#------------------------------------------------------------------------------#

setMethod("getRowsNames", signature(lp = "optObj_lpSolveAPI", i = "numeric"),

    function(lp, i) {

        rn <- dimnames(lp@oobj)[[1]][i]
        return(rn)

    }
)


#------------------------------------------------------------------------------#

setMethod("getColsNames", signature(lp = "optObj_lpSolveAPI", j = "numeric"),

    function(lp, j) {

        cn <- dimnames(lp@oobj)[[2]][j]
        return(cn)

    }
)


#------------------------------------------------------------------------------#


