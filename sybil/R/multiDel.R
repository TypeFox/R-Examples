#  multiDel.R
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
# Function: multDel
#
#
#
#

multiDel <- function(model, nProc = 2,
                     todo = "oneGeneDel",
                     del1 = NA, del2 = NA, ...) {
#                    #unlistResult = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    if (nProc < 2) {
        stop("argument nProc must be equal or greater than 2!")
    }

    #--------------------------------------------------------------#
    # split delX to a list of length numCo
    #--------------------------------------------------------------#
    spDel <- function(del) {
        if(is(del, "matrix")) {
            nd <- ncol(del)
            splitmat <- TRUE
        }
        else {
            nd <- length(del)
            splitmat <- FALSE
        }
        
        if (numCo > nd) {
            numCo <- nd
        }
        
        if (isTRUE(splitmat)) {
            gs  <- floor(seq(1, nd, length.out = numCo+1))
            dL <- vector(mode = "list", length = numCo)

            for (i in seq(along = gs[1:numCo])) {
                en <- ifelse(i == numCo, gs[i+1], gs[i+1]-1)
                dL[[i]] <- del[, gs[i]:en, drop = FALSE]
                #dL[[i]] <- c(gs[i], en)
            }
        }
        else {
            gs <- floor(seq(0, nd, length.out = numCo+1))
    
            # convert to factor
            spf <- cut(1:nd, gs)
            # use spf as factor for split()
            dL  <- split(del, spf) 
        }

        return(dL)
    }

    #--------------------------------------------------------------#


    # load library 'parallel'
    if (!requireNamespace("parallel", quietly = TRUE)) {
    	stop("package parallel not found.")
    }

#	unwanted conditioning for loading packages...
#    if(!isTRUE(require("parallel"))) {
#        stop("package parallel not found.")
#    }

    # number of cores
    ncore <- parallel::detectCores()
    
    numCo <- ifelse(nProc > ncore, as.integer(ncore), as.integer(nProc))


    #--------------------------------------------------------------#
    # split input into lists of size numCo

    if (any(is.na(del1))) {
        dL1 <- spDel(allGenes(model))
    }
    else {
        dL1 <- spDel(del1)
    }

    if (any(is.na(del2))) {
        dL2 <- as.list(rep(NA, length(dL1)))
    }
    else {
        if (length(del1) != length(del2)) {
            stop(paste("if argument del2 is not NA,",
                       "del1 and del2 must have same length!"))
        }
        dL2 <- spDel(del2)
        cdL <- vector("list", length(dL1))
        for (i in seq(along = cdL)) {
            cdL[[i]] <- c(dL1[i], dL2[i])
        }
    }


    #--------------------------------------------------------------#
    # run optimizations

    sol <- switch(todo,
        "oneGeneDel" = {
            parallel::mclapply(dL1,
                      function(x) oneGeneDel(model,
                                             geneList = x,
                                             verboseMode = 0, ...),
                      mc.cores = nProc)

        },
        "doubleGeneDel" = {
            parallel::mclapply(cdL, function(x) doubleGeneDel(model,
                                                     geneList1 = x[[1]],
                                                     geneList2 = x[[2]],
                                                     verboseMode = 0, ...),
                      mc.cores = nProc)
        },
        "oneFluxDel" = {
            parallel::mclapply(dL1, function(x) oneFluxDel(model,
                                                  react = x,
                                                  verboseMode = 0, ...),
                      mc.cores = nProc)
        },
        "doubleFluxDel" = {
            parallel::mclapply(cdL, function(x) doubleFluxDel(model,
                                                     react1 = x[[1]],
                                                     react2 = x[[2]],
                                                     verboseMode = 0, ...),
                      mc.cores = nProc)
        },
        "fluxVar" = {
            parallel::mclapply(dL1, function(x) fluxVar(model,
                                               react = x,
                                               verboseMode = 0, ...),
                      mc.cores = nProc)
        },
        "geneDeletion" = {
            parallel::mclapply(dL1, function(x) geneDeletion(model,
                                                    genes = x,
                                                    verboseMode = 0, ...),
                      mc.cores = nProc)
        },
        {
            stop("argument todo is not valid!")
        }
    )

#    if ( (isTRUE(unlistResult)) && (is(sol, "list")) ) {
#        ## maybe we can use something like unlist here
#        sv <- solver(sol[[1]])
#        nc <- lp_num_cols(sol[[1]])
#        nr <- lp_num_rows(sol[[1]])
#        of <- obj_function(sol[[1]])
#        np <- sum(mapply(num_of_prob, sol)) - (numCo-1)
#
#
#        # generate an index-vector of elements we want to use from sol
#        # every first element is the wild-type solution and this should be
#        # excluded, except for the first one
#
#        npL <- mapply(num_of_prob, sol, SIMPLIFY = FALSE, USE.NAMES = FALSE)
#        np  <- sum(unlist(npL)) - (numCo-1)
#
#        indL <- mapply(rep, npL, MoreArgs = list(x=TRUE), SIMPLIFY = FALSE)
#        
#        for (i in seq(along = indL)[-1]) {
#            indL[[i]][1] <- FALSE
#        }
#
#        newSol <- switch(todo,
#            "oneGeneDel" = {
#                optsol_geneDel(solver,
#                               nprob,
#                               lpdir,
#                               ncols,
#                               nrows,
#                               objf,
#                               fld,
#                               comb = 1)
#
#            },
#            "doubleGeneDel" = {
#                optsol_doublegenedel(solver,
#                                     nprob,
#                                     lpdir,
#                                     nrows,
#                                     ncols,
#                                     delrows,
#                                     delcols,
#                                     objf,
#                                     fld)
#            },
#            "oneFluxDel" = {
#                optsol_fluxdel(solver,
#                               nprob,
#                               lpdir,
#                               ncols,
#                               nrows,
#                               objf,
#                               fld,
#                               comb = 1)
#            },
#            "doubleFluxDel" = {
#                optsol_doublefluxdel(solver,
#                                     nprob,
#                                     lpdir,
#                                     nrows,
#                                     ncols,
#                                     delrows,
#                                     delcols,
#                                     objf,
#                                     fld)
#            },
#            "fluxVar" = {
#                optsol_fluxVar(solver,
#                               method,
#                               nprob,
#                               lpdir,
#                               ncols,
#                               nrows,
#                               objf,
#                               fld,
#                               rc)
#            },
#            "geneDeletion" = {
#                optsol_geneDel(solver,
#                               nprob,
#                               lpdir,
#                               ncols,
#                               nrows,
#                               objf,
#                               fld,
#                               comb = 1)
#            },
#            {
#                stop("argument todo is not valid!")
#            }
#        )
#    
#    }

    return(sol)

}
