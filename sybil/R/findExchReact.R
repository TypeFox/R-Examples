#  findExchReact.R
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
# Function: findExchReact
#
# The function findExchReact() is inspired by the function
# findExcRxns() contained in the COBRA Toolbox.
# The algorithm is the same.


findExchReact <- function(model) {

    if ( (!is(model, "modelorg")) &&
        !( (is(model, "Matrix")) || (is(model, "matrix")) ) ) {
        stop("needs an object of class modelorg, Matrix or matrix!")
    }

    if (is(model, "modelorg")) {
        St <- S(model)
    }
    else {
        St <- model
    }

    # columns with only one entry
    if(is(St, "Matrix")){
    	oneEntry <- colSums(St != 0)==1
    }
    else{
    	oneEntry <- apply(St, 2, function(x) sum(x != 0) == 1)
    }
    

    if (sum(oneEntry) > 0) {
        # exchange reactions -- with a -1 or 1
        
        if(is(St, "Matrix")){
			exchangeReact <- (colSums(St[ , oneEntry, drop = FALSE] == 1) == 1 | colSums(St[ , oneEntry, drop = FALSE]== -1) == 1)
		}
		else{
			exchangeReact <- apply(St[ , oneEntry, drop = FALSE], 2, function(x) (sum(x == 1) == 1) | (sum(x == -1) == 1))
		}
        
        # vector with the reaction id's of the exchange reactions
        ex <- c(1 : dim(St)[2])[oneEntry[exchangeReact]]

        # uptake reactions
        up <- NA
        if (is(model, "modelorg_irrev")) {
            up <- apply(St[,oneEntry], 2, function(x) (sum(x == 1) == 1))
        }

        if ((is(model, "modelorg")) && (!is(model, "modelorg_irrev"))) {
            up <- lowbnd(model)[ex] < 0
        }

        if (is(model, "modelorg")) {

            # get the row id's of S containing the non-zeros of the exchange reactions
            #exMet <- which(Matrix::rowSums(abs(S(model)[, ex])) == 1)

            # must be '> 0', because irreversible models (including exchange
            # reactions) have 2 entries per row, if the exchange reaction is
            # reversible
            exMet <- which(Matrix::rowSums(abs(S(model)[, ex, drop = FALSE])) > 0)

            # get the rows in the correct order
            # (if exchange reactions are not in main diagonal)

            ### as(S(model), "CsparseMatrix") ###
            metabolite <- exMet[S(model)[exMet, ex, drop = FALSE]@i+1]

            react <- new("reactId_Exch",
                         mod_id  = mod_id(model),                   
                         mod_key = mod_key(model),                   
                         rpnt    = ex,
                         rid     = react_id(model)[ex],
                         upt     = up,
                         mpnt    = metabolite,
                         mid     = met_id(model)[metabolite],
                         lb      = lowbnd(model)[ex],
                         ub      = uppbnd(model)[ex])
        }
        else {
            react <- oneEntry[exchangeReact]
        }
    }
    else {
        warning("no exchange reaction found")
        react <- NULL
    }

    return(react)

}
