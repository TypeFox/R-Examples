#  modelorg2ExPA.R
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
# Function: modelorg2ExPA
#
#
#
#
# original function buildExpaFileFromModelorg by C. Jonathan Fritzemeier


modelorg2ExPA <- function(model,
                          fname = NULL,
                          exIntReact = NULL,
                          filepath = ".",
                          suffix = "expa",
                          tol = SYBIL_SETTINGS("TOLERANCE")) {
	
    on.exit( closeAllConnections() )
    
    if (!is(model, "modelorg")) {
        stop("needs an object of class modelorg!")
    }
	
	#expa <- ifelse (missing(fname), mod_id(model), fname)
	if (is.null(fname)) {
	    expa <- paste(mod_id(model), suffix, sep = ".")
	}
	else {
	    expa <- fname
	}
	
    tofile <- file.path(filepath, expa)
    
    fh <- try(file(tofile, "wt"), silent = TRUE)
    
    if (is(fh, "try-error")) {
        warning("cannot write ExPA file!")
        fh <- FALSE
    }
    
    # exclude reactions
    if (is.null(exIntReact)) {
        exInt <- NA
    }
    else {
        exInt <- checkReactId(model, exIntReact)
    }
    
    # ------------------------------------------------------------------------ #
    # exchange reactions/excluded internal reactions
    # ------------------------------------------------------------------------ #

    exch <- findExchReact(model)
    ex <- logical(react_num(model))
    ex[react_pos(exch)] <- TRUE
    if (is(exInt, "reactId")) {
        ex[react_pos(exInt)] <- TRUE
    }

    # ------------------------------------------------------------------------ #
    # internal reactions
    # ------------------------------------------------------------------------ #
    
	cat("(Internal Fluxes)\n", file = fh)

	for (j in seq(along = react_id(model))) {
		
		if (isTRUE(ex[j])) {
		    next
		}
		
		# direction
		ri <- ifelse(isTRUE(react_rev(model)[j]), "R", "I")
		
		scol <- S(model)[ ,j]
		met  <- which(abs(scol) > tol)
        
        nint <- abs(scol[met] - round(scol[met])) > .Machine$double.eps^0.5
        if (any(nint)) {
            msg <- sprintf(ngettext(length(met_id(model)[met[nint]]),
                           "reaction %s contains %d non-integer stoichiometric coefficient, check metabolite %s",
                           "reaction %s contains %d non-integer stoichiometric coefficients, check metabolites\n\t%s"),
                           sQuote(react_id(model)[j]), length(met_id(model)[met[nint]]),
                           paste(sQuote(met_id(model)[met[nint]]), collapse = "\n\t"))
            warning(msg, call. = FALSE)
        }

		react <- character(length(met))
		
        # stoichiometric coefficient must not be rational, only integers
        # are valid for ExPA
		for (i in seq(along = met)) {

		    if (sign(scol[met[i]]) == 1) {
		        stoich <- paste("+", scol[met[i]], sep = "")
		    }
		    else {
		        stoich <- scol[met[i]]
		    }
		    
		    react[i] <- paste(stoich, met_id(model)[met[i]])
		}
		
		# reaction string
		rstr <- paste(react_id(model)[j],
		              ri,
		              paste(react, collapse = " "),
		              sep = "\t")

  		cat(rstr, "\n", sep = "", file = fh, append = TRUE)

	}
	

    # ------------------------------------------------------------------------ #
    # transport reactions
    # ------------------------------------------------------------------------ #

	cat("(Exchange Fluxes)\n", file = fh, append = TRUE)

    exInd <- react_pos(exch)
    
    for (i in seq(along = exInd)) {
    
        scol <- S(model)[ ,exInd[i]]
        met  <- which(scol != 0)
        
        if (length(met) != 1) {
            warning("error in reaction id", react_id(exch)[i])
            next
        }
        
        if ( (lowbnd(model)[exInd[i]] < 0) && (uppbnd(model)[exInd[i]] > 0) ) {
            exdir <- "Free"
        }
        else {
            if (lowbnd(model)[exInd[i]] < 0) {
                exdir <- "Input"
            }
            else {
                exdir <- "Output"
            }
        }
    
        estr <- paste(met_id(model)[met], exdir, sep = "\t")
        cat(estr, "\n", sep = "", file = fh, append = TRUE)
    
    }
	

    # ------------------------------------------------------------------------ #

    if (is(fh, "file")) {
        close(fh)
    }

    #--------------------------------------------------------------------------#
    # end
    #--------------------------------------------------------------------------#

    return(invisible(TRUE))

}

