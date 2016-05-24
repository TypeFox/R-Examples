#  modelorg2tsv.R
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
# Function: modelorg2tsv
#
# 
# 

modelorg2tsv <- function(model, prefix, suffix, extMetFlag = "b",
                         fielddelim = "\t", entrydelim = ", ",
                         makeClosedNetwork = FALSE,
                         onlyReactionList = FALSE,
                         minimalSet = FALSE,
                         fpath = SYBIL_SETTINGS("PATH_TO_MODEL"),
                         ...) {

    ## on.exit( closeAllConnections() )
    
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

    # validate model structure before writing
    validObject(model)

    # filenames
    if (missing(prefix)) {
        prefix <- gsub("\\s+", "_", mod_id(model))
    }
    
    if (missing(suffix)) {
        suffix <- switch(fielddelim,
            "\t" = { "tsv" },
            ";"  = { "csv" },
            ","  = { "csv" },
            "|"  = { "dsv" },
                   { "dsv" }
        )
    }
    
    fnameR <- paste(paste(prefix, "react", sep = "_"), suffix, sep = ".")
    fnameM <- paste(paste(prefix, "met",   sep = "_"), suffix, sep = ".")
    fnameD <- paste(paste(prefix, "desc",  sep = "_"), suffix, sep = ".")

    # path to output file
    tsvfileR <- file.path(fpath, fnameR)
    tsvfileM <- file.path(fpath, fnameM)
    tsvfileD <- file.path(fpath, fnameD)


    #--------------------------------------------------------------------------#
    # reactions list
    #--------------------------------------------------------------------------#

    # create reaction strings
    rstr <- .createReactionString(model,
                                          makeClosedNetwork,
                                          entrydelim,
                                          extMetFlag)

    if (isTRUE(onlyReactionList)) {
        write.table(x = data.frame(equation = rstr$equat),
                    row.names = FALSE, file = tsvfileR, sep = fielddelim, ...)
    }
    else if (isTRUE(minimalSet)) {
        write.table(x = data.frame(
                    abbreviation = react_id(model),
                    equation     = rstr$equat,
                    lowbnd       = lowbnd(model),
                    uppbnd       = uppbnd(model),
                    obj_coef     = obj_coef(model)),
            row.names = FALSE, file = tsvfileR, sep = fielddelim, ...)
    }
    else {

        # some optional entries
        if (length(gpr(model)) != react_num(model)) {
            gp <- character(react_num(model))
        }
        else {
            gp <- gpr(model)
        }
        
        if (nrow(subSys(model)) != react_num(model)) {
            sys <- character(react_num(model))
        }
        else {
            sys_tmp    <- subSys(model)
            ssys_names <- colnames(sys_tmp)

            if ( (length(ssys_names) == 1) && is.na(ssys_names) ) {
                sys <- character(react_num(model))
            }
            else {
                sys <- apply(sys_tmp, 1,
                             function(x) paste(ssys_names[x], collapse = ", "))
            }
        }

        write.table(x = data.frame(
                    abbreviation = react_id(model),
                    name         = react_name(model),
                    equation     = rstr$equat,
                    reversible   = rstr$revers,
                    compartment  = rstr$compat,
                    lowbnd       = lowbnd(model),
                    uppbnd       = uppbnd(model),
                    obj_coef     = obj_coef(model),
                    rule         = gp,
                    subsystem    = sys),
            row.names = FALSE, file = tsvfileR, sep = fielddelim, ...)
    }


    #--------------------------------------------------------------------------#
    # metabolites list
    #--------------------------------------------------------------------------#

    if ( (!isTRUE(onlyReactionList)) && (!isTRUE(minimalSet)) ) {

        metunq <- sort(unique(rstr$metab))
        
        mpos <- lapply(metunq, function(x) rstr$metab %in% x)
        
        # metabolite names
        metNames <- lapply(mpos, function(x) unique(met_name(model)[x]))
        metNames <- unlist(lapply(metNames, paste, collapse = entrydelim))
        
        # metabolite compartments
        metCompart <- lapply(mpos,
                             function(x) mod_compart(model)[met_comp(model)[x]])
        metCompart <- unlist(lapply(metCompart, paste, collapse = entrydelim))
        
        write.table(x = data.frame(
                    abbreviation = metunq,
                    name         = metNames,
                    compartment  = metCompart),
            row.names = FALSE, file = tsvfileM, sep = fielddelim, ...)

    }
    

    #--------------------------------------------------------------------------#
    # model description
    #--------------------------------------------------------------------------#

    if ( (!isTRUE(onlyReactionList)) && (!isTRUE(minimalSet)) ) {

        # get id's of metabolites in different compartments
        # (one per compartment)
        metDiffComp <- match(mod_compart(model),
                             mod_compart(model)[met_comp(model)])
        
        metAbbrevComp <- character(length(metDiffComp))

        # get the compartment abbreviations
        metALl <- grepl("^.+(\\[\\w+\\])$", met_id(model)[metDiffComp])
        
        metAbbrevComp[metALl]  <- sub("^.+(\\[\\w+\\])$", "\\1",
                                      met_id(model)[metDiffComp[metALl]])

        metAbbrevComp[!metALl] <- mod_compart(model)[!metALl]
        
        # generate output format
        ma <- paste(metAbbrevComp, collapse = entrydelim)
        mc <- paste(mod_compart(model), collapse = entrydelim)

        write.table(x = data.frame(
                    name         = mod_name(model),
                    id           = mod_id(model),
                    description  = mod_desc(model),
                    compartment  = mc,
                    abbreviation = ma,
                    Nmetabolites = met_num(model),
                    Nreactions   = react_num(model),
                    Ngenes       = length(allGenes(model)),
                    Nnnz         = Snnz(model)),
            row.names = FALSE, file = tsvfileD, sep = fielddelim, ...)

    }


    #--------------------------------------------------------------------------#
    # end
    #--------------------------------------------------------------------------#

    return(TRUE)

}
