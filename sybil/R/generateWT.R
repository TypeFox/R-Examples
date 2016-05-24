#  generateWT.R
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
# Function: .generateWT
#
#
#

.generateWT <- function(model, react = NULL, lb = NULL, ub = NULL, ...) {

#    ca <- match.call()
#    
#    if ("solver" %in% names(ca)) {

#        # It is necessary to test, whether argument solver has a variable,
#        # or not. If it is a variable, it needs to be evaluated.
#        testslv <- tryCatch(eval(parse(text = ca["solver"])),
#                            error = function(e) e)
#        if (is(testslv, "simpleError")) {
#            slv <- as.character(ca["solver"])
#        }
#        else {
#            slv <- as.character(testslv)
#        }
#        
#    }
#    else {
#        slv <- SYBIL_SETTINGS("SOLVER")
#    }
	
	# setting ... parameters into list
	ca <- list(...)
	if(is.null(ca$solver)){
		ca$solver <- SYBIL_SETTINGS("SOLVER")
	}

    me <- checkDefaultMethod(solver = ca$solver,
                             method = "NA",
                             probType = "lp",
                             loadPackage = FALSE)
    
    ca$solver <- me[["sol"]]
    ca$method <- me[["met"]]
    ca$solverParm <- as.data.frame(NA)
    
    ca$object <- model
    if(is.null(ca$algorithm)) ca$algorithm <- "fba"
    ca$react <- react
    ca$ub <- ub
    ca$lb <- lb

    if (is(react, "list")) {
        message("calculating fba solutions ... ", appendLF = FALSE)
        suppressMessages({
        	ca$lpdir <- rep("max", length(react))
        	ca$verboseMode <- 0
        	
        	tmp <- do.call(optimizer, ca)
#            tmp <- optimizer(model, algorithm = algorithm,
#                             lpdir = rep("max", length(react)),
#                             react = react, lb = lb, ub = ub, verboseMode = 0, 
#                             solver = me[["sol"]], method = me[["met"]],
#                             solverParm = as.data.frame(NA), ...)
        })
        message("OK")
    }
    else {
    	ca$retOptSol <- FALSE
    	ca$lpdir <- "max"
    	tmp <- do.call("optimizeProb", ca)
#        tmp <- optimizeProb(model,
#                            react = react, lb = lb, ub = ub,
#                            retOptSol = FALSE,
#                            algorithm = algorithm,
#                            lpdir = "max",
#                            solver = me[["sol"]],
#                            method = me[["met"]],
#                            solverParm = as.data.frame(NA), ...)
    }

    return(tmp)

}
