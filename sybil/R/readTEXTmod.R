#  readTEXTmod.R
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
# Function: readTEXTmod
#
# 
# 

readTEXTmod <- function(filename,
                        description,
                        def_bnd = SYBIL_SETTINGS("MAXIMUM")) {

if ( file.exists(filename) == FALSE ) {
    stop( c("cannot open file ", filename) )
}

if (missing(description)) {
    description <- filename
}


#------------------------------------------------------------------------------#
#                              some functions                                  #
#------------------------------------------------------------------------------#

parse_metabolites <- function(reactands, fwbw) {

    new_met <- character(0)

    stoichnum <- gregexpr("\\b[0-9]+(\\.[0-9]*)?\\b", reactands)

    # the metabolites -> \w{2,} so we do not get the "+"
    #print(ed_prod[1])
    metabolites <- gregexpr("([0-9]+(\\.[0-9]*)?\\s)?\\w{2,}\\b", reactands)
    #print(metabolites)


    for (l in 1:length(metabolites[[1]])) {
        metabolite <- substr(reactands, metabolites[[1]][l], metabolites[[1]][l] + attr(metabolites[[1]], "match.length")[l] - 1)
        print(metabolite)
        
        
        stoich_met <- unlist(strsplit(metabolite, " ", fixed = TRUE))
        if (length(stoich_met) == 1) {
            stoich <- 1 * fwbw
            met <- stoich_met
        }
        else {
            stoich <- as.numeric(stoich_met[1]) * fwbw
            met <- stoich_met[2]
        }

        # exclude external metabolites
        if (regexpr("ex$", met) != -1) {
            next
        }

        met_number <- which(met_id %in% met)
        print(met_number)
        # new metabolite
        if (length(met_number) == 0) {
            # add a new row to S
            if (length(met_id) > 0) {
                S <- rbind(S, numeric(dim(S)[2]))
            }
            met_id <- c(met_id, met)
            S[length(met_id), num_react] <- stoich
        }
        
        # existing metabolite
        if (length(met_number) == 1) {
            S[met_number, num_react] <- stoich
        }

        # existing metabolite
        if (length(met_number) > 1) {
            stop(paste("Error in line ", i, ": check metabolite ", met, "!", sep = ""))
        }

    }

}

#------------------------------------------------------------------------------#

parse_genes <- function(gene) {

    if (length(gene) == 0) {
        return(list(gpr_str = "", rules_str = "", gene_vec = ""))
    }

    # for the gpr slot and gprRules slot
    gene_pos <- which(allGenes %in% gene)
    if (length(gene) > 1) {
        gpr_string <- paste("(", paste(gene, sep = "", collapse = " and "), ")")
        rules_string <- paste("x[", gene_pos, "]", sep = "")
        rules_string <- paste("(", paste(rules_string, sep = "", collapse = " & "), ")")
    }
    else {
        gpr_string <- gene
        rules_string <- paste("x[", gene_pos, "]", sep = "")
    }
#     print(gpr_string)
#     print(rules_string)
    return(list(gpr_str = gpr_string, rules_str = rules_string, gene_vec = gene))

}

#------------------------------------------------------------------------------#

get_stoich_met <- function(reactand) {

    stoich_met <- unlist(strsplit(reactand, " ", fixed = TRUE))

    if (length(stoich_met) == 2) {
        if(is.na(as.numeric(stoich_met[1]))) {
        print(reactand)
            stop(paste("Something went wrong in line", i, "with a stoichiometric number!"))
        }
        
        stoich <- as.numeric(stoich_met[1])
        met <- stoich_met[2]
    }
    if (length(stoich_met) == 1) {
        stoich <- 1
        met <- stoich_met[1]
    }
    if (length(stoich_met) > 2) {
        stop(paste("Something went wrong in line", i))
    }

    return(list(stoich = stoich, met = met))

}


#------------------------------------------------------------------------------#
#                            reading the model                                 #
#------------------------------------------------------------------------------#

message("reading the model ...", appendLF = FALSE)

ModelTmp <- readLines(filename)
num_lines <- length(ModelTmp)

# The new model
Model <- modelorg(filename, filename)
mod_desc(Model) <- filename

message(" OK.")


# needed variables

met_id <- character(0)
react_id <- character(0)
react_name <- character(0)
react_rev <- logical(0)
num_react <- 0
num_met <- 0

ext_met <- character(0)

S <- matrix(0)

genes <- list()
rules <- character(0)
gpr <- character(0)
allGenes <- character(0)

lowbnd <- numeric(0)
uppbnd <- numeric(0)

external_metabolites <- character(0)

message("parsing reactions ...", appendLF = FALSE)

last_reactions <- character(0)
reactions <- strsplit(ModelTmp, "\t|[ ]{4}")
#reactions <- strsplit(ModelTmp, "\t")

for (i in 1 : num_lines) {

# print(i)
# print(length(gpr))
# print(length(genes))
# print("")

    current_genes <- character(0)

    # ignore empty lines
    if (length(reactions[[i]]) == 0) {next}

    # Lines beginning with a # are treated as comment lines.
    if (substr(reactions[[i]][1], 0, 1) == "#") {next}

    # Do here somethin with the first line
    
    
    #print(reactions[[i]][3])
    #print("last")
    #print(last_reactions)


    # reactions we already have
    # We should do here the same as in doubleGenes(): check S for that task.
    # That will be more safer to typos.
    lrids <- which(last_reactions %in% reactions[[i]][3])
    #print(lrids)
    if (length(lrids) == 1) {

        # react_name slot
        tmp_name <- unlist(strsplit(react_name[lrids], ":", fixed = TRUE))
        tmp_name[1] <- paste(tmp_name[1], reactions[[i]][2], sep = "|")
        react_name[lrids] <- paste(tmp_name[1], tmp_name[2], sep = ":")  


        current_genes <- unlist(strsplit(reactions[[i]][1], "/", fixed = TRUE))
    
        # allGenes slot
        allGenes <- unlist(union(allGenes, current_genes))

        gene_list <- parse_genes(current_genes)

        # genes slot
        genes[[lrids]] <- c(genes[[lrids]], gene_list$gene_vec)
 
        if (gpr[lrids] == "") {
            gpr[lrids] <- gene_list$gpr_str
            rules[lrids] <- gene_list$rules_str
        }
        else {
            gpr[lrids] <- paste(gpr[lrids], gene_list$gpr_str, sep = " or ")
            rules[lrids] <- paste(rules[lrids], gene_list$rules_str, sep = " | ")
        }

        next
    }

    if (length(lrids) > 1) {
        stop(paste("Error!", i))
    }

    num_react <- num_react + 1

    exchange_reaction <- FALSE

    # add a new column to S
    if (dim(S)[1] != 1) {
        S <- cbind(S, numeric(dim(S)[1]))
    }

    # react_id slot
    react_id <- c(react_id, paste("v", num_react, sep = ""))    

    # react_name slot
    react_name <- c(react_name, paste(reactions[[i]][2], reactions[[i]][3], sep = ": "))    

    current_genes <- unlist(strsplit(reactions[[i]][1], "/", fixed = TRUE))

    # allGenes slot
    allGenes <- unlist(union(allGenes, current_genes))


    gene_list <- parse_genes(current_genes)

    # genes slot
    genes[[length(genes) + 1]] <- gene_list$gene_vec
    gpr <- c(gpr, gene_list$gpr_str)
    rules <- c(rules, gene_list$rules_str)

    # parse reactions

    # get the reaction arrow (maybe there is a more elegant way to do this)
    #dir_coord <- regexpr("[<>-]+", reactions[[i]][3])
    dir_coord <- regexpr("(<-+)|(-+>)|(<-+>)", reactions[[i]][3])
    direction <- substr(reactions[[i]][3], dir_coord[1], dir_coord[1] + attr(dir_coord, "match.length") - 1)

#    print(reactions[[i]][3])

    ed_prod <- unlist(strsplit(reactions[[i]][3], direction, fixed = TRUE))
    if (length(ed_prod) != 2) {
        #print(reactions[[i]])
        #print(ed_prod)
        stop(paste("Reaction in line", i, "is strange!"))
    }
#    print(ed_prod)


    # check direction of reaction
    check <- logical(2)

    # backward
    if (regexpr("^<", direction) != -1) {
        check[1] <- TRUE
    }

    # forward
    if (regexpr(">$", direction) != -1) {
        check[2] <- TRUE
    }

    # reaction is reversible
    if (all(check == TRUE)) {
        react_rev <- c(react_rev, TRUE)
    }
    else {
        if (all(check == FALSE)) {
            stop(paste("Maybe no reaction arrow in line", i, "?"))
        }
        
        # backward
        if (check[1] == TRUE) {
            react_rev <- c(react_rev, FALSE)
            ed_prod <- rev(ed_prod)
        }
        
        # forward
        if (check[2] == TRUE) {
            react_rev <- c(react_rev, FALSE)
        }
    }

    ## Here some work is still to do --> create functions

    # The Educts

    
    reactands <- unlist(strsplit(ed_prod[1], " + ", fixed = TRUE))
    reactands <- sub("^\\s|\\s+$", "", reactands)

    for (l in 1:length(reactands)) {

        st_me_list <- get_stoich_met(reactands[l])

        stoich <- st_me_list$stoich * -1
        met <- st_me_list$met

        if (regexpr("xt$", met) != -1) {
            ext_met <- c(ext_met, met)
        }

        met_number <- which(met_id %in% met)
        #print(met_number)
        # new metabolite
        if (length(met_number) == 0) {
            num_met <- num_met + 1
            # add a new row to S
            if (length(met_id) > 0) {
                S <- rbind(S, numeric(dim(S)[2]))
            }
            met_id <- c(met_id, met)
            S[num_met, num_react] <- stoich
        }
        
        # existing metabolite
        if (length(met_number) == 1) {
            S[met_number, num_react] <- stoich
        }

        # existing metabolite
        if (length(met_number) > 1) {
            stop(paste("Error in line ", i, ": check metabolite ", met, "!", sep = ""))
        }


    }


    # The Products

    reactands <- unlist(strsplit(ed_prod[2], " + ", fixed = TRUE))
    reactands <- sub("^\\s|\\s+$", "", reactands)

    for (l in 1:length(reactands)) {

        st_me_list <- get_stoich_met(reactands[l])

        stoich <- st_me_list$stoich
        met <- st_me_list$met
        
        if (regexpr("xt$", met) != -1) {
            ext_met <- c(ext_met, met)
        }

        met_number <- which(met_id %in% met)
        #print(met_number)
        # new metabolite
        if (length(met_number) == 0) {
            num_met <- num_met + 1
            # add a new row to S
            if (length(met_id) > 0) {
                S <- rbind(S, numeric(dim(S)[2]))
            }
            met_id <- c(met_id, met)
            S[num_met, num_react] <- stoich
        }
        
        # existing metabolite
        if (length(met_number) == 1) {
            S[met_number, num_react] <- stoich
        }

        # existing metabolite
        if (length(met_number) > 1) {
            stop(paste("Error in line ", i, ": check metabolite ", met, "!", sep = ""))
        }


    }

    # the lower and upper bounds
    
    # lower

    if (length(reactions[[i]]) > 3) {
        lowbnd <- c(lowbnd, as.numeric(reactions[[i]][4]))
    }
    else {
        if (all(check == TRUE)) {
            lowbnd <- c(lowbnd, def_bnd * -1)
        }
        else {
            lowbnd <- c(lowbnd, 0)
        }
    }


    # upper
    
    if (length(reactions[[i]]) > 4) {
        uppbnd <- c(uppbnd, as.numeric(reactions[[i]][5]))
    }
    else {
        uppbnd <- c(uppbnd, def_bnd)
    }

    last_reactions <- c(last_reactions, reactions[[i]][3])

}

message(" OK.")


# exchange reactions

if (length(ext_met) > 0) {

    message("exchange reactions ...", appendLF = FALSE)
    
    ext_met <- unique(ext_met)
    num_ext_met <- length(ext_met)
    ext_met_ids <- which(met_id %in% ext_met)
    
    
    S <- cbind(S, matrix(0, num_met, num_ext_met))
    #apply(as.matrix(c(1:num_ext_met)), 1, function(x) S[ext_met_ids[x], (num_react + x)] <- -1)
    for (i in 1:num_ext_met) {
        S[ext_met_ids[i], (num_react + i)] <- -1
    }
    
    num_react <- num_react + num_ext_met
    
    react_id <- c(react_id, paste("b", c(1:num_ext_met), sep = ""))
    react_name <- c(react_name, paste(ext_met, "exchange", sep = "_"))
    react_rev <- c(react_rev, rep(TRUE, num_ext_met))
    
    uppbnd <- c(uppbnd, rep(def_bnd, num_ext_met))
    lowbnd <- c(lowbnd, rep(0, num_ext_met))

    gpr <- c(gpr, rep("", num_ext_met))
    rules <- c(rules, rep("", num_ext_met))
    genes[length(genes):num_react] <- ""
    
    
    message(" Ok.")

}



# rxnGeneMat

message("GPR mapping ...", appendLF = FALSE)
genes_ind <- lapply(genes, function(x) allGenes %in% x)

rxnGeneMat <- matrix(0, num_react, length(allGenes))

for (i in seq(along = genes_ind)) {
    rxnGeneMat[i, genes_ind[[i]]] <- 1
}
message(" Ok.")




met_id(Model)     <- met_id
react_id(Model)   <- react_id
react_name(Model) <- react_name
react_rev(Model)  <- react_rev
react_num(Model)  <- as.integer(num_react)
met_num(Model)    <- as.integer(num_met)

S(Model)          <- S

genes(Model)      <- genes
gprRules(Model)   <- rules
gpr(Model)        <- gpr
allGenes(Model)   <- allGenes
rxnGeneMat(Model) <- rxnGeneMat

obj_coef(Model)   <- integer(num_react)

uppbnd(Model)     <- uppbnd
lowbnd(Model)     <- lowbnd
# uppbnd(Model) <- rep(def_bnd, num_react)
# lowbnd(Model) <- integer(num_react)
# lowbnd(Model)[react_rev] <- def_bnd * -1

return(Model)



}
