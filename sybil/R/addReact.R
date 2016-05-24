#  addReact.R
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
# Function: addReact
#
#
# The function addReact() is inspired by the function
# addReaction() contained in the COBRA Toolbox.
# The algorithm is (more or less) the same.


addReact <- function(model,
                     id,
                     met,
                     Scoef,
                     reversible = FALSE,
                     lb = 0,
                     ub = SYBIL_SETTINGS("MAXIMUM"),
                     obj = 0,
                     subSystem = NA,
                     gprAssoc = NA,
                     reactName = NA,
                     metName = NA,
                     metComp = NA) {

  
    # ------------------------------------------------------------------------ #
    # check arguments
    # ------------------------------------------------------------------------ #

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    if (length(met) != length(Scoef)) {
        stop("arguments 'met' and 'Scoef' must have the same length")
    }

    if (length(id) > 1) {
        stop("add/change one reaction")
    }

    if ( ( (ub > 0) && (lb < 0) ) && (!isTRUE(reversible)) ) {
        Crev <- TRUE
        warning(paste("'lb' and 'ub' are signed different,",
                      "set reversible to 'TRUE'"))
    }
    else {
       Crev <- reversible
    }


    # ------------------------------------------------------------------------ #
    # check, if we need to add columns and/or rows
    # ------------------------------------------------------------------------ #

    # reaction
    colInd <- match(id, react_id(model))
    addCol <- FALSE
    nCols  <- react_num(model)
    if (is.na(colInd)) {
        # new reaction
        colInd <- react_num(model) + 1
        addCol <- TRUE
        nCols  <- nCols + 1
    }


    # metabolites
    rowInd <- match(met, met_id(model))
    
    newM     <- which(is.na(rowInd))
    nRows    <- met_num(model)          # number of rows in the model
    nNewRows <- length(newM)            # number of new rows
    addRow   <- FALSE
    
    for (i in seq(along = newM)) {
        addRow   <- TRUE
        nRows    <- nRows + 1
        rowInd[newM[i]] <- nRows
    }
    

    if ( (isTRUE(addCol)) || (isTRUE(addRow)) ) {    

        # -------------------------------------------------------------------- #
        # make a new model
        # -------------------------------------------------------------------- #

        # -------------------------------------------------------------------- #
        # data structures

        newmet_num      <- met_num(model)
        newmet_id       <- met_id(model)
        newmet_name     <- met_name(model)
        newmet_comp     <- met_comp(model)
        newmet_single   <- met_single(model)
        newmet_de       <- met_de(model)

        newreact_num    <- react_num(model)
        newreact_rev    <- react_rev(model)
        newreact_id     <- react_id(model)
        newreact_name   <- react_name(model)
        newreact_single <- react_single(model)
        newreact_de     <- react_de(model)
        newlowbnd       <- lowbnd(model)
        newuppbnd       <- uppbnd(model)
        newobj_coef     <- obj_coef(model)

        newgprRules     <- gprRules(model)
        newgenes        <- genes(model)
        newgpr          <- gpr(model)
        newallGenes     <- allGenes(model)
        newrxnGeneMat   <- rxnGeneMat(model)
        newsubSys       <- subSys(model)

        newS            <- S(model)

    
        if (isTRUE(addRow)) {
            
            # new number of metabolites
            newmet_num  <- nRows
            
            # new metabolite id's
            newmet_id   <- append(met_id(model), met[newM])

            # new metabolite names
            if (any(is.na(metName))) {
                newmet_name <- append(met_name(model), met[newM])
            }
            else {
                newmet_name <- append(met_name(model), metName[newM])
            }

            # new metabolite compartments
            if (any(is.na(metComp))) {
                newmet_comp <- append(met_comp(model), rep(NA, nNewRows))
            }
            else {
                if (is(metComp, "numeric")) {
                    newmet_comp <- append(met_comp(model), metComp[newM])
                }
                else {
                    newmet_comp <- append(met_comp(model),
                                          match(metComp[newM],
                                                mod_compart(model)))
                }
            }

            # singleton and dead end metabolites (not checked!)
            newmet_single <- append(met_single(model), rep(NA, nNewRows))
            newmet_de     <- append(met_de(model),     rep(NA, nNewRows))

            # new rows in stoichiometric matrix
            newRows <- Matrix::Matrix(0,
                                      nrow = nNewRows,
                                      ncol = react_num(model))
            newS <- Matrix::rBind(newS, newRows)
        }
    
        if (isTRUE(addCol)) {                        # we add at most one column
            # new number of reactions
            newreact_num  <- nCols
            
            # new reaction id
            newreact_id   <- append(react_id(model), id)

            # new reaction name
            if (is.na(reactName)) {
                newreact_name <- append(react_name(model), id)
            }
            else {
                newreact_name <- append(react_name(model), reactName)
            }

            # reaction contains singleton or dead end metabolites (not checked!)
            newreact_single <- append(react_single(model), NA)
            newreact_de     <- append(react_de(model),     NA)

            # reversibility, lower and upper bounds, objective coefficient
            newreact_rev <- append(react_rev(model), Crev)
            newlowbnd    <- append(lowbnd(model),    lb)
            newuppbnd    <- append(uppbnd(model),    ub)
            newobj_coef  <- append(obj_coef(model),  obj)

            # new column in stoichiometric matrix
            newS <- cBind(newS, rep(0, nrow(newS)))
            
            # subsystems
            if (any(is.na(subSystem))) {
            	ss <- subSys(model)
            	if(ncol(ss)==0){ # if no subSys defined, rbind (see else) failed
            		dim(ss) <- c(nrow(ss)+1, ncol(ss))
            		newsubSys <- ss
            	}
            	else {
            		newsubSys <- rBind(ss, rep(FALSE, ncol(subSys(model))))
            	}
            }
            else {
                if (is(subSystem, "logical")) {
                    newsubSys <- rBind(subSys(model), subSystem)
                }
                else {
                    nSubsRow  <- colnames(subSys(model)) %in% subSystem
                    newsubSys <- rBind(subSys(model), nSubsRow)
                }
            }


            # gpr association
            if (ncol(rxnGeneMat(model)) > 0) {
                newrxnGeneMat   <- rBind(rxnGeneMat(model),
                                         rep(FALSE, ncol(rxnGeneMat(model))))
            }
            else { #if (nrow(rxnGeneMat(model)) > 0) {
                newrxnGeneMat <- rxnGeneMat(model)
                dim(newrxnGeneMat) <- c(nrow(newrxnGeneMat)+1,
                                        ncol(newrxnGeneMat))
            }
            # do above else always.

            if ( (is.na(gprAssoc)) || (gprAssoc == "") ) {
                if ((length(gprRules(model)) > 0)) {
                    newgprRules     <- append(gprRules(model), "")
                    newgenes        <- append(genes(model), "")
                    newgpr          <- append(gpr(model), "")
                }
            }
            else {
                gene_rule <- .parseBoolean(gprAssoc)
            
                geneInd <- match(gene_rule$gene, allGenes(model))
            
                # indices of new genes
                new_gene <- which(is.na(geneInd))

                # if we have new gene(s), add a column in rxnGeneMat and
                # gene name(s) to allGenes
                if (length(new_gene) > 0) {
                    newallGenes <- append(allGenes(model),
                                          gene_rule[["gene"]][new_gene])

                    # update geneInd
                    geneInd <- match(gene_rule[["gene"]], newallGenes)

                    # if we have an empty modelorg object, we need to
                    # initialize rxnGeneMat
                    if (ncol(newrxnGeneMat) == 0) {
                        newrxnGeneMat <- Matrix::Matrix(FALSE,
                                                        nCols, max(geneInd))
                    }
                    else {
                        for (i in seq(along = gene_rule[["gene"]][new_gene])) {
	    					newrxnGeneMat <- cBind(newrxnGeneMat,
		    								   rep(FALSE, nrow(newrxnGeneMat)))
			    		}
					}
                }

                # rxnGeneMat
                newrxnGeneMat[nCols, geneInd] <- TRUE
 
                # new rule
                newgpr <- append(gpr(model), gprAssoc)

                # genes per reaction
                newgenes <- append(genes(model), list(gene_rule$gene))
                newrule  <- gene_rule$rule

                for (j in 1 : length(geneInd)) {
                    pat  <- paste("x(", j, ")", sep = "")
                    repl <- paste("x[", geneInd[j], "]", sep = "")
    
                    newrule <- gsub(pat, repl, newrule, fixed = TRUE)
                }

                newgprRules <- append(gprRules(model), newrule)
            }
        }
        
        # values for stoichiometric matrix
        newS[ , colInd]      <- 0    
        newS[rowInd, colInd] <- Scoef    
        
#        for (i in seq(along = rowInd)) {
#            newS[rowInd[i], colInd] <- Scoef[i]
#        }
    
        # -------------------------------------------------------------------- #
        # new model
        # -------------------------------------------------------------------- #
        
        if (is(model, "modelorg_irrev")) {
            mod_out <- modelorg_irrev(mod_id(model), mod_name(model))
            irrev(mod_out)     <- TRUE
            matchrev(mod_out)  <- append(matchrev(model), 0)
            
            revReactId <- max(irrev2rev(model))+1
            irrev2rev(mod_out) <- append(irrev2rev(model), revReactId)
            rev2irrev(mod_out) <- rbind(rev2irrev(model), c(nCols, nCols))
        }
        else {
            mod_out <- modelorg(mod_id(model), mod_name(model))
        }
        
        mod_desc(mod_out)    <- mod_desc(model)
        mod_compart(mod_out) <- mod_compart(model)


        met_num(mod_out)      <- as.integer(newmet_num)
        met_id(mod_out)       <- newmet_id
        met_name(mod_out)     <- newmet_name
        met_comp(mod_out)     <- as.integer(newmet_comp)
        met_single(mod_out)   <- newmet_single
        met_de(mod_out)       <- newmet_de

        react_num(mod_out)    <- as.integer(newreact_num)
        react_rev(mod_out)    <- newreact_rev
        react_id(mod_out)     <- newreact_id
        react_name(mod_out)   <- newreact_name
        react_single(mod_out) <- newreact_single
        react_de(mod_out)     <- newreact_de
        lowbnd(mod_out)       <- newlowbnd
        uppbnd(mod_out)       <- newuppbnd
        obj_coef(mod_out)     <- newobj_coef

        gprRules(mod_out)     <- newgprRules
        genes(mod_out)        <- newgenes
        gpr(mod_out)          <- newgpr
        allGenes(mod_out)     <- newallGenes
        rxnGeneMat(mod_out)   <- newrxnGeneMat
        subSys(mod_out)       <- newsubSys

        S(mod_out)            <- newS

    }
    else {
    
        # -------------------------------------------------------------------- #
        # modify old model
        # -------------------------------------------------------------------- #
        
        mod_out <- model
        
        react_rev(mod_out)[colInd] <- Crev
        lowbnd(mod_out)[colInd]    <- lb
        uppbnd(mod_out)[colInd]    <- ub
        obj_coef(mod_out)[colInd]  <- obj
        S(mod_out)[ , colInd]      <- 0
        S(mod_out)[rowInd, colInd] <- Scoef

    }
    
    
    check <- validObject(mod_out, test = TRUE)

    if (check != TRUE) {
        msg <- paste("Validity check failed:", check, sep = "\n    ")
        warning(msg)
    }

    return(mod_out)

}

