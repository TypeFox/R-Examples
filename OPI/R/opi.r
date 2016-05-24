#
# Choose which OPI implementation to run and distribute calls accordingly.
#
# This would all have been nicer in an OO style, with each implementation
# being a subclass of an opi class, but I don't think it can be in R.
# The OPI standard doesn't want users employing exactly the same function 
# no matter what the underlying implementation, and so there cannot be 
# extra parameters to create method signatures for different classes.
# Similarly, some implementations use exactly the same method signatures,
# again which will confuse R, I think. Anyway, if I am wrong, sorry about that.
# What I've done (use a list of implementations and then use a global
# integer to index them) works and makes sense to the non-OO person.
#
# Author: Andrew Turpin    (aturpin@unimelb.edu.au)
# Date: June 2012
# Modified:    Sep 2014 - added Octopus 600
#           16 Dec 2014 - added Kowa AP 7000
#
# Copyright 2012 Andrew Turpin and Jonathan Denniss
# This program is part of the OPI (http://perimetry.org/OPI).
# OPI is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

.OpiEnv <- new.env(size=1)

################################################################################
# A list of available OPI implementations for chooseOpi to choose from, and 
# the opi* functions to index using .OpiEnv$chooser.
################################################################################
opi.implementations <- list(
    list(
        name="KowaAP7000",
        opiInitialize    = "kowaAP7000.opiInitialize",
        opiClose         = "kowaAP7000.opiClose",
        opiSetBackground = "kowaAP7000.opiSetBackground",
        opiQueryDevice   = "kowaAP7000.opiQueryDevice",
        opiPresent       = "kowaAP7000.opiPresent"
    ),
    list(
        name="Octopus900",
        opiInitialize    = "octo900.opiInitialize",
        opiClose         = "octo900.opiClose",
        opiSetBackground = "octo900.opiSetBackground",
        opiQueryDevice   = "octo900.opiQueryDevice",
        opiPresent       = "octo900.opiPresent"
    ),
    list(
        name="Octopus900F310",
        opiInitialize    = "octo900.opiInitialize",
        opiClose         = "octo900.opiClose",
        opiSetBackground = "octo900.opiSetBackground",
        opiQueryDevice   = "octo900.opiQueryDevice",
        opiPresent       = "octo900.opiPresentF310"
    ),
    list(
      name="Octopus600",
      opiInitialize    = "octo600.opiInitialize",
      opiClose         = "octo600.opiClose",
      opiSetBackground = "octo600.opiSetBackground",
      opiQueryDevice   = "octo600.opiQueryDevice",
      opiPresent       = "octo600.opiPresent"
    ),
    list(
        name="SimNo",
        opiInitialize    = "simNo.opiInitialize",
        opiClose         = "simNo.opiClose",
        opiSetBackground = "simNo.opiSetBackground",
        opiQueryDevice   = "simNo.opiQueryDevice",
        opiPresent       = "simNo.opiPresent"
    ),
    list(
        name="SimYes",
        opiInitialize    = "simYes.opiInitialize",
        opiClose         = "simYes.opiClose",
        opiSetBackground = "simYes.opiSetBackground",
        opiQueryDevice   = "simYes.opiQueryDevice",
        opiPresent       = "simYes.opiPresent"
    ),
    list(
        name="SimHenson",
        opiInitialize    = "simH.opiInitialize",
        opiClose         = "simH.opiClose",
        opiSetBackground = "simH.opiSetBackground",
        opiQueryDevice   = "simH.opiQueryDevice",
        opiPresent       = "simH.opiPresent"
    ),
    list(
        name="SimGaussian",
        opiInitialize    = "simG.opiInitialize",
        opiClose         = "simG.opiClose",
        opiSetBackground = "simG.opiSetBackground",
        opiQueryDevice   = "simG.opiQueryDevice",
        opiPresent       = "simG.opiPresent"
    ),
    list(
        name="SimHensonRT",
        opiInitialize    = "simH_RT.opiInitialize",
        opiClose         = "simH_RT.opiClose",
        opiSetBackground = "simH_RT.opiSetBackground",
        opiQueryDevice   = "simH_RT.opiQueryDevice",
        opiPresent       = "simH_RT.opiPresent"
    )
)

################################################################################
# Input parameters
#   opiImplementation  Either "Octopus900", "HEP", "SimHenson", "SimGaussian"
#                      If NULL, prints a list of possible values. Returns TRUE.
# Side effect
#   Sets .OpiEnv$chooser
#
# Returns
#   TRUE     If successful
#   FALSE    Otherwise
################################################################################
chooseOpi <- function(opiImplementation) { 
    possible <- unlist(lapply(opi.implementations, "[", "name"))

        #
        # If NULL, print the list of possible
        #
    if (missing(opiImplementation))
        opiImplementation <- NULL
    if (is.null(opiImplementation)) {
        print(possible)
        return(TRUE)
    }

        #
        # Warn about the one unimplemented one
        #
    if (opiImplementation == "HEP") {
        # require(rHEP)
        warning("Have not implemented chooseOPI(HEP)")
        return(FALSE)
    } 

        #
        # Check OPIOctopus900 package exists
        #
    #if ((opiImplementation == "Octopus900") && !require(OPIOctopus900)) {
    #    cat("***********************************************************************\n")
    #    cat("* You cannot choose the Octopus900 OPI without installing the package *\n")
    #    cat("* OPIOctopus900, which is available with permission from HAAG-STREIT. *\n")
    #    cat("***********************************************************************\n")
    #    stop("Get the Octopus900 package")
    #}

        #
        # Find the index in opi.implementations
        #
    i <- which(opiImplementation == possible)
    if (length(i) == 0) {
        assign("chooser", NA, envir=.OpiEnv)
        warning(paste("chooseOpi() cannot find opiImplementation",opiImplementation))
        return(FALSE)
    } else {
        assign("chooser", i, envir=.OpiEnv)
        return(TRUE)
    }
}


####################################################################################
# Simply send the opi*() call to the right implementation
####################################################################################
opiDistributor <- function(method, ...) {
    if (!exists("chooser", where=.OpiEnv) || is.na(.OpiEnv$chooser)) {
        msg <- "You have not chosen a valid OPI implementation. Use chooseOpi()"
        warning(msg)
        return(msg)
    }
    toCall <- opi.implementations[[.OpiEnv$chooser]][[method]]
    allowedArgs <- names(formals(toCall))
    haveArgs    <- names(list(...))
#print(paste("Allowed args: ", allowedArgs))
#print(paste("Have args: ", haveArgs))
    argsToPass  <- intersect(allowedArgs, haveArgs)
    argsNotPassed  <- setdiff(haveArgs, argsToPass)

    if (length(argsNotPassed) > 0)
        warning(paste(method, "Ignored argument ", argsNotPassed, "\n"))
#print(paste("Passing args: ", argsToPass))
    do.call(toCall, list(...)[argsToPass])
}

opiPresent        <- function(stim,nextStim=NULL,...) { opiDistributor("opiPresent", stim=stim, nextStim=nextStim, ...) }

opiInitialize     <- function(...) { opiDistributor("opiInitialize", ...) }
opiInitialise     <- function(...) { opiDistributor("opiInitialize", ...) }

opiSetBackground  <- function(...) { opiDistributor("opiSetBackground", ...) }

opiQueryDevice    <- function(...) { opiDistributor("opiQueryDevice", ...) }

opiClose          <- function(...) { opiDistributor("opiClose", ...) }
