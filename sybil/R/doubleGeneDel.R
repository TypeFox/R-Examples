#  doubleGeneDel.R
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
# Function: doubleGeneDel
#
# This function performs a "double gene deletion analysis".
# In each iteration one gene is switched of (vi = 0)
# and the objective function will be computed.
#
# The function doubleGeneDel() is inspired by the function
# doubleGeneDeletion() contained in the COBRA Toolbox.


doubleGeneDel <- function(model, geneList1, geneList2, lb = NULL, ub = NULL,
                          allComb = FALSE, exLethal = TRUE,
                          tol = SYBIL_SETTINGS("TOLERANCE"),
                          checkOptSolObj = FALSE, ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }


    # maybe we need here something like checkGeneId

    if (missing(geneList1)) {
        geneList1 <- allGenes(model)
    }
    else {
        if (is.na(all(match(geneList1, allGenes(model))))) {
            stop("check genelist1!")
        }
    }

    if (missing(geneList2)) {
        geneList2 <- allGenes(model)
    }
    else {
        if (is.na(all(match(geneList2, allGenes(model))))) {
            stop("check genelist2!")
        }
    }

    if (length(geneList1) < 1) {
        stop("Argument 'geneList1' must contain at least one gene!")
    }
    if (length(geneList2) < 1) {
        stop("Argument 'geneList2' must contain at least one gene!")
    }

    #geneList1 <- sort(match(geneList1, allGenes(model)))
    #geneList2 <- sort(match(geneList2, allGenes(model)))
    geneList1 <- match(geneList1, allGenes(model))
    geneList2 <- match(geneList2, allGenes(model))


    num_geneList1 <- length(geneList1)
    num_geneList2 <- length(geneList2)

    if ((!isTRUE(allComb)) && (num_geneList1 != num_geneList2)) {
        stop("geneList1 and geneList2 must have the same length")
    }

    num_genes <- length(allGenes(model))


#------------------------------------------------------------------------------#
#                           remove lethal genes                                #
#------------------------------------------------------------------------------#

    if (isTRUE(exLethal)) {
    
        # all different genes from list1 and list2
        unGenes <- sort(unique(c(geneList1, geneList2)))

        ca <- match.call()
        if ("solver" %in% names(ca)) {
            # It is necessary to test, whether argument solver has a variable,
            # or not. If it is a variable, it needs to be evaluated.
            testslv <- tryCatch(eval(parse(text = ca["solver"])),
                                error = function(e) e)
            if (is(testslv, "simpleError")) {
                slv <- as.character(ca["solver"])
            }
            else {
                slv <- as.character(testslv)
            }
        }
        else {
            slv <- SYBIL_SETTINGS("SOLVER")
        }

        wtobj <- optimizeProb(model, retOptSol = FALSE,
                              algorithm = "fba", solver = slv, lpdir = "max")

        unSol <- oneGeneDel(model, geneList = unGenes, solver = slv,
                            lpdir = rep("max", length(unGenes)),
                            fld = "none", algorithm = "fba")

        # solution id of lethal genes
        letid <- lethal(unSol, wtobj$obj)

        # check wether some of the 'essential' genes result in
        # an unsuccessfull solution
        statNok <- checkStat(unSol[letid])

        lethal <- letid
        #lethal <- letid[-statNok]

        # lethalGeneIds contain the ids of the lethal genes
        lethalGeneIds <- unGenes[lethal]

        l1 <- which(geneList1 %in% lethalGeneIds)
        l2 <- which(geneList2 %in% lethalGeneIds)

        if (isTRUE(allComb)) {
        
            if (length(l1) > 0) {
                geneList1 <- geneList1[-l1]
            }
            if (length(l2) > 0) {
                geneList2 <- geneList2[-l2]
            }    

        }
        else {
            l12 <- sort(unique(c(l1, l2)))
            if (length(l12) > 0) {
                geneList1 <- geneList1[-l12]
                geneList2 <- geneList2[-l12]
            }
        }
        remove(unSol)
     }
     else {
         lethalGeneIds <- as.integer(NA)
     }


#------------------------------------------------------------------------------#
#               calculate the number of optimizations (num_opt)                #
#------------------------------------------------------------------------------#

    # m <- outer(mystring, mystring, paste, sep="")
    # m[upper.tri(m)]

    
    # use merge() or expand.grid() or combn() here!!

    if (isTRUE(allComb)) {
  
        # Compute Boolean matrix with TRUE in the upper triangonal
        # (the maximum number of comparisons)
        tmpMAT <- upper.tri(matrix(nrow = num_genes,
                                   ncol = num_genes), diag = FALSE)

        # The next step is, to compute the differences between the two fields.
        # If an element of react1 is not in react2, we have to set the
        # corresponding row to TRUE and vice versa.
        tmpDIFF <- setdiff(geneList2, geneList1)
        #print(tmpDIFF)
    
        if (length(tmpDIFF) != 0) {
            tmpMAT[,tmpDIFF] <- TRUE
        }
    
        tmpDIFF <- setdiff(geneList1, geneList2)
        #print(tmpDIFF)
    
        if (length(tmpDIFF) != 0) {
            tmpMAT[tmpDIFF,] <- TRUE
        }

        # Exclude everything from tmpMAT we do not need.
        #tmpMAT <- tmpMAT[geneList1, geneList2, drop = FALSE]
      
    }
    else {
  
        tmpMAT <- matrix(FALSE, nrow = num_genes, ncol = num_genes)
        #diag(tmpMAT) <- TRUE
        for (i in seq(along = geneList1)) {
            tmpMAT[geneList1[i], geneList2[i]] <- TRUE
        }
    }

    geneList11 <- geneList1
    geneList22 <- geneList2
    
    geneList1 <- unique(geneList1)
    geneList2 <- unique(geneList2)

    # Exclude everything from tmpMAT we do not need.
    tmpMAT <- tmpMAT[geneList1, geneList2, drop = FALSE]
  
    # The number of TRUE's in tmpMAT is equal to the number of optimizations
    num_opt <- sum(tmpMAT == TRUE)


    rownames(tmpMAT) <- geneList1
    colnames(tmpMAT) <- geneList2
    #print(tmpMAT)
    #print(num_opt)  


#------------------------------------------------------------------------------#
#                               run optimization                               #
#------------------------------------------------------------------------------#

    deletions <- which(tmpMAT == TRUE, arr.ind = TRUE)
 
    kogenesID <- cbind(geneList1[deletions[,"row"]],
                       geneList2[deletions[,"col"]])
    kogenes   <- lapply(seq_len(nrow(kogenesID)), function(x) kogenesID[x, ])

    fd <- .generateFluxdels(model, kogenes)

    if (is.null(lb)) {
        lb <- rep(0, length(kogenes))
    }
    else {
        if (length(lb) != length(kogenes)) {
            stop("lb must be of length ", length(kogenes))
        }
    }
    if (is.null(ub)) {
        ub <- rep(0, length(kogenes))
    }
    else {
        if (length(ub) != length(kogenes)) {
            stop("ub must be of length ", length(kogenes))
        }
    }

    sol <- optimizer(model = model,
                     react = fd[["react"]],
                     lb    = lb,
                     ub    = ub,
                     ...)


    # ------------------------------------------------------------------------ #

    optsol <- new("optsol_genedel")
    opt <- makeOptsolMO(model, sol)
    as(optsol, "optsol_optimizeProb") <- opt
    
    chlb(optsol)      <- as.numeric(lb)
    chub(optsol)      <- as.numeric(ub)
    dels(optsol)      <- matrix(allGenes(model)[kogenesID], ncol = 2)
    fluxdels(optsol)  <- fd[["fd"]]
    hasEffect(optsol) <- fd[["heff"]]

    if (isTRUE(checkOptSolObj)) {
        checkOptSol(optsol, onlywarn = TRUE)
    }

    return(optsol)

}

