#  optimizer.R
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


################################################
# Function: optimizer
#
#
#


optimizer <- function(model, react,
                      lb = NULL,
                      ub = NULL,
                      obj_coef = NULL,
                      lpdir = NULL,
                      algorithm = SYBIL_SETTINGS("ALGORITHM"),
                      mtfobj = NULL,
                      setToZero = FALSE,
                      rebuildModel = FALSE,
                      fld = "none",
                      prCmd = NA, poCmd = NA,
                      prDIR = NULL, poDIR = NULL,
                      verboseMode = 2,
                      ...) {


    stopifnot(length(fld) == 1)

    #--------------------------------------------------------------------------#
    # verboseMode

    on.exit(expr = {
        if (exists("logObj")) {
            logClose(logObj) <- NA
        }
    } )

    # start logging
    logObj <- sybilLog(filename = "",
                       loglevel = -1,
                       verblevel = verboseMode)


    #--------------------------------------------------------------------------#
    # check arguments

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg")
    }

    if (!is(react, "list")) {
        stop("needs an object of class list")
    }

    # number of optimizations
    nObj  <- length(react)

    # parameters modifying the model
    if (!is.null(lb)) {
        stopifnot((length(lb) == nObj || nrow(lb) == nObj),
                  (is(lb, "numeric")  || is(lb, "list") || is(lb, "matrix")))
    }
    if (!is.null(ub)) {
        stopifnot((length(ub) == nObj || nrow(ub) == nObj),
                  (is(ub, "numeric")  || is(ub, "list") || is(ub, "matrix")))
    }
    if (!is.null(obj_coef)) {
        stopifnot((length(obj_coef) == nObj || nrow(obj_coef) == nObj),
                  (is(obj_coef, "numeric") ||
                   is(obj_coef, "list")    ||
                   is(obj_coef, "matrix")))
    }
    if (!is.null(lpdir)) {
        stopifnot(length(lpdir) == nObj, is(lpdir, "character"))
    }
    
    # flux distribution
    if (isTRUE(fld)) {
        fdist <- "all"
    }
    else if (identical(fld, FALSE)) {
        fdist <- "none"
    }
    else {
        fdist <- fld
    }


    #--------------------------------------------------------------------------#
    # convenient function

    gEl <- function(el, num, pos) {

        if (is.null(el)) {
            return(el)
        }
        
        stopifnot(length(pos) == 1)

        if (is.list(el)) {
            ret <- el[[pos]]
        }
        else if (is.matrix(el)) {
            ret <- el[pos, , drop = TRUE]
        }
        else {
            ret <- rep(el[pos], num)
        }
    
        stopifnot(length(ret) == num)
        return(ret)
    
    } 
#    gEl <- function(el, num) {
#    
#        if (is.null(el)) {
#            return(el)
#        }
#        
#        stopifnot(length(el) == 1)
#
#        if (is.list(el)) {
#            ret <- unlist(el)
#            stopifnot(length(el) == num)
#        }
#        else {
#            ret <- rep(el, num)
#        }
#    
#        return(ret)
#    
#    } 


    #--------------------------------------------------------------------------#
    # prepare problem object
    #--------------------------------------------------------------------------#

    if (algorithm == "mtf") {
        if (fdist == "none") {
            fdist <- "fluxes" 
        }
        if (is.null(mtfobj)) {
            lpmod  <- sysBiolAlg(model, algorithm = "mtf",
                                 react = react, lb = lb, ub = ub, ...)
        }
        else {
            stopifnot(is(mtfobj, "numeric"), length(mtfobj) == nObj)
            lpmod  <- sysBiolAlg(model, algorithm = "mtf", wtobj = mtfobj, ...)
        }
    }
    else {
        lpmod  <- sysBiolAlg(model, algorithm = algorithm, ...)
    }

    # check, if we use an algorithm performing genetic perturbations
    pert <- checkAlgorithm(algorithm, "pert")


    #--------------------------------------------------------------------------#
    # data structures for simulation results
    #--------------------------------------------------------------------------#

    obj   <- numeric(nObj)
    mobj  <- numeric(nObj)
    ok    <- integer(nObj)
    stat  <- integer(nObj)
    flux <- switch(fdist,
        "all" = {
            Matrix::Matrix(0, nrow = nc(lpmod), ncol = nObj)
        },
        "fluxes" = {
            Matrix::Matrix(0, nrow = length(fldind(lpmod)), ncol = nObj)
        },
        {
            NA
        }
    )


    #--------------------------------------------------------------------------#
    # pre and post processing

    runPrPl  <- logical(nObj)
    runPoPl  <- logical(nObj)
    runPrPcn <- 1
    runPoPcn <- 1

    if (all(!is.na(prCmd))) {
        do_pr  <- TRUE
        prPcmd <- NULL
        runPrP <- .doInRound(prDIR, nObj)
        prPpa  <- vector(mode = "list", length = length(runPrP))
        runPrPl[runPrP] <- TRUE
    }
    else {
        do_pr <- FALSE
    }
    if (all(!is.na(poCmd))) {
        do_po  <- TRUE
        poPcmd <- NULL
        runPoP <- .doInRound(poDIR, nObj)
        poPpa  <- vector(mode = "list", length = length(runPoP))
        runPoPl[runPoP] <- TRUE
    }
    else {
        do_po <- FALSE
    }


#------------------------------------------------------------------------------#
#                             optimizations                                    #
#------------------------------------------------------------------------------#

    message("calculating ", nObj, " optimizations ... ", appendLF = FALSE)
    if (verboseMode > 1) { cat("\n") }
    if (verboseMode == 2) {
        progr <- .progressBar()
        #progr <- txtProgressBar(min = 2, max = nObj, initial = 2, style = 3)
    }

    logOptimizationTH(logObj)


    fi <- fldind(lpmod)
    objcTMP <- integer(react_num(model))


    for (i in 1:nObj) {

        if (verboseMode == 2) {
            progr <- .progressBar(i, nObj, progr)
            #setTxtProgressBar(progr, i)
        }

        # pre/post processing
        if (isTRUE(runPrPl[i])) {
            prCmd_tmp <- prCmd
            did_pr    <- TRUE
        }
        else {
            prCmd_tmp <- NA
            did_pr    <- FALSE
        }

        if (isTRUE(runPoPl[i])) {
            poCmd_tmp <- poCmd
            did_po    <- TRUE
        }
        else {
            poCmd_tmp <- NA
            did_po    <- FALSE
        }

        # solution i
        if (isTRUE(rebuildModel)) {
            sol <- optimizeProb(model,
                                react = react[[i]],
#                                lb = rep(lb[i], length(react[[i]])),
#                                ub = rep(ub[i], length(react[[i]])),
#                                obj_coef = obj_coef[i],
#                                lb = gEl(lb[i], length(react[[i]])),
#                                ub = gEl(ub[i], length(react[[i]])),
#                                obj_coef = gEl(obj_coef[i], length(react[[i]])),
                                lb = gEl(lb, length(react[[i]]), i),
                                ub = gEl(ub, length(react[[i]]), i),
                                obj_coef = gEl(obj_coef, length(react[[i]]), i),
                                lpdir = lpdir[i],
                                #lpdir = getObjDir(problem(lpmod)),
                                retOptSol = FALSE,
                                prCmd = prCmd_tmp, poCmd = poCmd_tmp,
                                prCil = runPrPcn, poCil = runPoPcn,
                                algorithm = algorithm, ...)
        }
        else {
            if (algorithm == "mtf") {
                changeMaxObj(lpmod, i)
            }

            sol <- optimizeProb(lpmod,
                                react = react[[i]],
#                                lb = rep(lb[i], length(react[[i]])),
#                                ub = rep(ub[i], length(react[[i]])),
#                                obj_coef = obj_coef[i],
#                                lb = gEl(lb[i], length(react[[i]])),
#                                ub = gEl(ub[i], length(react[[i]])),
#                                obj_coef = gEl(obj_coef[i], length(react[[i]])),
                                lb = gEl(lb, length(react[[i]]), i),
                                ub = gEl(ub, length(react[[i]]), i),
                                obj_coef = gEl(obj_coef, length(react[[i]]), i),
                                lpdir = lpdir[i],
                                prCmd = prCmd_tmp, poCmd = poCmd_tmp,
                                prCil = runPrPcn, poCil = runPoPcn)
        }

        #obj[i]  <- sum(obj_coef(model) * sol$fluxes[fi])
        #mobj[i] <- crossprod(obj_coef(model), sol$fluxes[fi])
        #obj[i]  <- sol$obj
        ok[i]   <- sol$ok
        stat[i] <- sol$stat
        if (fdist == "none") {
            if (isTRUE(pert)) {
                #obj[i] <- crossprod(obj_coef(model), sol$fluxes[fi])
                if (is.null(obj_coef)) {
                    obj[i] <- crossprod(obj_coef(model), sol$fluxes[fi])
                }
                else {
                    objcTMP[react[[i]]] <- gEl(obj_coef, length(react[[i]]), i)
                    obj[i] <- crossprod(objcTMP, sol$fluxes[fi])
                }
            }
            else {
                obj[i] <- sol$obj
            }
        }
        else {
            obj[i] <- sol$obj
            if (fdist == "fluxes") {
                flux[,i] <- sol$fluxes[fi]
            }
            else {
                flux[,i] <- sol$fluxes
            }
        }

        # pre/post processing
        if (isTRUE(did_pr)) {
            if ( (runPrPcn == 1) && (is.null(prPcmd)) ) {
                prPcmd  <- cmd(sol$preP)
            }
            prPpa[[runPrPcn]] <- pa(sol$preP)
            runPrPcn <- runPrPcn+1
            did_pr <- FALSE
        }
        if (isTRUE(did_po)) {
            if ( (runPoPcn == 1) && (is.null(poPcmd)) ) {
                poPcmd  <- cmd(sol$postP)
            }
            poPpa[[runPoPcn]] <- pa(sol$postP)
            runPoPcn <- runPoPcn+1
            did_po <- FALSE
        }

        logOptimization(logObj, sol$ok, sol$stat, obj[i], lpdir[i], obj_coef[i], react[[i]], i)

        remove(sol)
        #close(progr)
    }

    message("OK")


#------------------------------------------------------------------------------#
#                               save the results                               #
#------------------------------------------------------------------------------#

    # slot fldind
    if (fdist == "fluxes") {
        fli <- 1:length(fi)
    }
    else if(fdist == "none") {
        fli <- NA
    }
    else {
        fli <- fi
    }

    # pre and post processing
    if (isTRUE(do_pr)) {
        prAna <- ppProc(prPcmd)
        pa(prAna) <- prPpa
        ind(prAna) <- runPrP
    }
    else {
        prAna <- NULL
    }

    if (isTRUE(do_po)) {
        poAna <- ppProc(poPcmd)
        pa(poAna) <- poPpa
        ind(poAna) <- runPoP
    }
    else {
        poAna <- NULL
    }

    # solution list
    optsol <- list(solver       = solver(problem(lpmod)),
                   method       = method(problem(lpmod)),
                   algorithm    = algorithm(lpmod),
                   lp_num_cols  = nc(lpmod),
                   lp_num_rows  = nr(lpmod),
                   obj          = obj,
                   ok           = ok,
                   stat         = stat,
                   lp_dir       = factor(getObjDir(problem(lpmod))),
                   fldind       = fli,
                   fluxdist     = fluxDistribution(flux),
                   prAna        = prAna,
                   poAna        = poAna,
                   alg_par      = alg_par(lpmod))


#------------------------------------------------------------------------------#

    if (isTRUE(setToZero)) {
        do_again <- checkSolStat(stat, solver(problem(lpmod)))
        num_new  <- length(do_again)
        optsol[["obj"]][do_again] <- as.numeric(0)

        message("setting ", num_new, " objective values to zero")

        for (i in seq(along = do_again)) {
            logOptimization(logObj,
                            optsol[["ok"]][do_again[i]],
                            optsol[["stat"]][do_again[i]], 0,
                            lpdir[[do_again[i]]], obj_coef[[do_again[i]]],
                            react[[do_again[i]]], do_again[i])
        }
    }


#------------------------------------------------------------------------------#
#                           return solution object                             #
#------------------------------------------------------------------------------#

    delProb(problem(lpmod))
    remove(lpmod)

    logFoot(logObj)  <- TRUE
    logClose(logObj) <- NA

    return(optsol)

}



