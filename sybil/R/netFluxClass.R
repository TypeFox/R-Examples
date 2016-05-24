#  netFluxClass.R
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


# netFluxClass


#------------------------------------------------------------------------------#
#                            class definitions                                 #
#------------------------------------------------------------------------------#

setClass("netFlux",
    representation(
        uptake       = "logical",
        product      = "logical",
        unused       = "logical",
        react_id     = "character",
        rate         = "numeric"
    ),
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

getNetFlux <- function(rates, tol = SYBIL_SETTINGS("TOLERANCE")) {

    id  <- names(rates)
    names(rates) <- NULL

    upt <- rates < tol * -1
    prd <- rates > tol
    uusd <- logical(length(rates))
    uusd[!upt] <- TRUE
    uusd[!prd] <- TRUE

    nf <- new("netFlux",
              react_id = as.character(id),
              rate = as.numeric(rates),
              uptake = upt, product = prd, unused = uusd)

    return(nf)
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# react_id
setMethod("react_id", signature(object = "netFlux"),
          function(object) {
              return(object@react_id)
          }
)

setReplaceMethod("react_id", signature = (object = "netFlux"),
                 function(object, value) {
                     object@react_id <- value
                     return(object)
                 }
)


# rate
setMethod("rate", signature(object = "netFlux"),
          function(object) {
              return(object@rate)
          }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#



# consider using sprintf here
setMethod("show", signature(object = "netFlux"),
    function(object) {
        ri <- react_id(object)
        rr <- rate(object)
        cat("uptake reaction rates (absolute values):\n")
        cat(sprintf(" %-20s%10f\n", ri[object@uptake], abs(rr[object@uptake])), sep = "")
        cat("\nexcretion reaction rates (absolute values):\n")
        cat(sprintf(" %-20s%10f\n", ri[object@product], rr[object@product]), sep = "")
        cat("\nunused exchange reactions [abs(rate) < ", SYBIL_SETTINGS("TOLERANCE"), "]:\n", sep = "")
        print(ri[object@unused])
    }
)


# length of an object of class optsol
setMethod("length", signature = signature(x = "netFlux"),
          function(x) {
              return(length(rate(x)))
          }
)

