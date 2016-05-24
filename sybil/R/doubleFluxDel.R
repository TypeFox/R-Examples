#  doubleFluxDel.R
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
# Function: doubleFluxDel
#
# This function performs a "double gene deletion analysis".
# In each iteration one gene is switched of (vi = 0)
# and the objective function will be computed.


doubleFluxDel <- function(model, react1, react2, lb = NULL, ub = NULL,
                          allComb = FALSE, exex = FALSE,
                          checkOptSolObj = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    ex <- 0

    if (missing(react1)) {
        ex <- 1
        react1 <- checkReactId(model, 1:react_num(model))
    }
	else {
        if (!is(react1, "reactId")) {
            react1 <- checkReactId(model, react1)
        }
    }

    if (missing(react2)) {
        ex <- 2
        react2 <- checkReactId(model, 1:react_num(model))
    }
    else {
        if (!is(react2, "reactId")) {
            react2 <- checkReactId(model, react2)
        }
    }

    react1 <- (react_pos(react1)) # removed sort here
    react2 <- (react_pos(react2)) # removed sort here


#------------------------------------------------------------------------------#
#                          remove exchange reactions                           #
#------------------------------------------------------------------------------#

    if ((ex == 2) && (exex == TRUE)) {

        exch <- findExchReact(model)
        exch <- react_pos(exch)

        react1 <- react1[-exch]
        react2 <- react2[-exch]

    }

    num_react1 <- length(react1)
    num_react2 <- length(react2)
 
    if ((!isTRUE(allComb)) && (num_react1 != num_react2)) {
        stop("react1 and react2 must have the same length")
    }


#------------------------------------------------------------------------------#
#               calculate the number of optimizations (num_opt)                #
#------------------------------------------------------------------------------#

#     print(react1)
#     print(react2)
  
    if (isTRUE(allComb)) {
		# if allComb is used, duplicated entries in react1 or react2 dont make 
		# sense:
		react1 <- unique(react1)
		react2 <- unique(react2)
		
		
        # Compute Boolean matrix with TRUE in the upper triangonal
        # (the maximum number of comparisons)
        tmpMAT <- upper.tri(matrix(nrow = react_num(model),
                                   ncol = react_num(model)))
    
        # The next step is, to compute the differences between the two fields.
        # If an element of react1 is not in react2, we have to set the
        # corresponding row to TRUE and vice versa.
        tmpDIFF <- setdiff(react2, react1)
        # print(tmpDIFF)
    
        if (length(tmpDIFF) != 0) {
            tmpMAT[,tmpDIFF] <- TRUE
        }
    
        tmpDIFF <- setdiff(react1, react2)
        # print(tmpDIFF)
    
        if (length(tmpDIFF) != 0) {
            tmpMAT[tmpDIFF,] <- TRUE
        }
    
        # Exclude everything from tmpMAT we do not need.
        tmpMAT <- tmpMAT[react1, react2]
      
        # We need this, if the resulting matrix is no longer 2 dimensional
        if (is.null(dim(tmpMAT))) {
            tmpMAT <- as.matrix(tmpMAT, byrow = TRUE)
            if (num_react1 < num_react2) {
                tmpMAT <- t(tmpMAT)
            }
        }
        # The number of TRUE's in tmpMAT is equal to the number of optimizations
        num_opt <- sum(tmpMAT == TRUE)
        
        
        rownames(tmpMAT) <- react1
	    colnames(tmpMAT) <- react2
        
		deletions <- which(tmpMAT == TRUE, arr.ind = TRUE)
		koreactID <- cbind(react1[deletions[,"row"]],
				  		   react2[deletions[,"col"]])

    }
    else {
  		 koreactID <- cbind(react1, react2)
#        tmpMAT <- matrix(FALSE, nrow = react_num(model),
#                                ncol = react_num(model))
#        for (i in 1:num_react1) {
#            tmpMAT[react1[i], react2[i]] <- TRUE
#        }
#        #tmpMAT[geneList1, geneList2] <- TRUE
#        browser()
#        num_opt <- num_react1
#        tmpMAT <- tmpMAT[react1, react2]
    }

    koreact   <- lapply(seq_len(nrow(koreactID)), function(x) koreactID[x, ])

    # The number of TRUE's in tmpMAT is equal to the number of optimizations
    # print(num_opt)

#------------------------------------------------------------------------------#
#                               run optimization                               #
#------------------------------------------------------------------------------#
    
    if (is.null(lb)) {
        lb <- rep(0, length(koreact))
    }
    else {
        if (length(lb) != length(koreact)) {
            stop("lb must be of length ", length(koreact))
        }
    }
    if (is.null(ub)) {
        ub <- rep(0, length(koreact))
    }
    else {
        if (length(ub) != length(koreact)) {
            stop("ub must be of length ", length(koreact))
        }
    }

    sol <- optimizer(model = model, lb = lb, ub = ub, react = koreact, ...)


    # ------------------------------------------------------------------------ #

    optsol <- new("optsol_fluxdel")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt
    
    chlb(optsol) <- as.numeric(lb)
    chub(optsol) <- as.numeric(ub)
    dels(optsol) <- matrix(react_id(model)[koreactID], ncol = 2)

    if (isTRUE(checkOptSolObj)) {
        checkOptSol(optsol, onlywarn = TRUE)
    }

    return(optsol)

}

