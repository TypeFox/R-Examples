#  modelorg2text.R
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
# Function: modelorg2text
#
# 
# 

modelorg2text <- function(model, prefix, suffix, extMetFlag = "b",
                         fielddelim = "\t", genedelim = "/",
                         makeClosedNetwork = FALSE,
                         fpath = SYBIL_SETTINGS("PATH_TO_MODEL"),
                         ...) {

    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }

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
    
    fname <- paste(prefix, suffix, sep = ".")

    # path to output file
    textfile <- file.path(fpath, fname)


    #--------------------------------------------------------------------------#
    # some functions
    #--------------------------------------------------------------------------#

    prepareRuleStrings <- function(rule, bp, no) {
    
    
    }


    #--------------------------------------------------------------------------#
    # reactions list
    #--------------------------------------------------------------------------#

    rstr <- .createReactionString(model, makeClosedNetwork)

    # rule    abbreviation    equation    lowbnd    uppbnd    obj_coef
    
    gprFRule <- gsub("&&?|and",    "AND", gpr(model), perl = TRUE)
    gprFRule <- gsub("\\|\\|?|or", "OR",  gprFRule,   perl = TRUE)
    
    gpr_list <- strsplit(gprFRule, "", fixed = TRUE)
    
    check_bp <- mapply(.check_brackets, gpr_list, SIMPLIFY = TRUE)
   
    if ( sum(check_bp) != length(check_bp) ) {
       warning(paste("Wrong gpr rules detected, setting to \"\". ",
                     "Check rule(s) no. ",
                     paste((1:react_num(model))[!check_bp], collapse = ", "),
                     ".", sep = ""))
       gpr_list[!check_bp] <- ""
    }
    
    gpr_bp  <- mapply(.bracket_pairs, gpr_list, SIMPLIFY = FALSE)
    
    gpr_str <- mapply(prepareRuleStrings,
                      gprFRule, gpr_bp, c(1:react_num(model)),
                      SIMPLIFY = FALSE)
    
    
#    write.table(x = data.frame(
#                    equation = rstr$equat,
#                    lowbnd       = lowbnd(model),
#                    uppbnd       = uppbnd(model),
#                    obj_coef     = obj_coef(model),
#                    ),
#        row.names = FALSE, file = tsvfileR, sep = fielddelim, ...)



    #--------------------------------------------------------------------------#
    # end
    #--------------------------------------------------------------------------#

    return(TRUE)

}
