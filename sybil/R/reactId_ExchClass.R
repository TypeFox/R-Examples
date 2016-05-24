#  reactId_ExchClass.R
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


# reactId_ExchClass


#------------------------------------------------------------------------------#
#                   definition of the class reactId_Exch                       #
#------------------------------------------------------------------------------#

setClass("reactId_Exch",
         representation(
              uptake  = "logical",
              met_pos = "integer",
              met_id  = "character",
              lowbnd  = "numeric",
              uppbnd  = "numeric"
         ),
         contains = "reactId",
         validity = .validreactId_Exch
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "reactId_Exch",
          definition = function(.Object, mod_id, mod_key = "",
                                rpnt, rid, upt, mpnt, mid, lb, ub) {

              if ( (!missing(mod_id))  ||
                   (!missing(mod_key)) ||
                   (!missing(rpnt))    ||
                   (!missing(rid))     ||
                   (!missing(upt))     ||
                   (!missing(mpnt))    ||
                   (!missing(mid))     ||
                   (!missing(lb))      ||
                   (!missing(ub)) ) {
              
                  tmp <- new("reactId",
                             mod_id  = mod_id,
                             mod_key = mod_key,
                             pnt     = rpnt,
                             id      = rid)
                  as(.Object, "reactId_Exch") <- tmp
              
                  .Object@uptake  <- as.logical(upt)
                  .Object@met_pos <- as.integer(mpnt)
                  .Object@met_id  <- as.character(mid)
                  .Object@lowbnd  <- as.numeric(lb)
                  .Object@uppbnd  <- as.numeric(ub)
                  validObject(.Object)
              }

              return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# position
setMethod("met_pos", signature(object = "reactId_Exch"),
          function(object) {
              return(object@met_pos)
          }
)

setReplaceMethod("met_pos", signature = (object = "reactId_Exch"),
                 function(object, value) {
                     object@met_pos <- value
                     return(object)
                 }
)


# met_id
setMethod("met_id", signature(object = "reactId_Exch"),
          function(object) {
              return(object@met_id)
          }
)

setReplaceMethod("met_id", signature = (object = "reactId_Exch"),
                 function(object, value) {
                     object@met_id <- value
                     return(object)
                 }
)


# lowbnd
setMethod("lowbnd", signature(object = "reactId_Exch"),
          function(object) {
              return(object@lowbnd)
          }
)

setReplaceMethod("lowbnd", signature = (object = "reactId_Exch"),
                 function(object, value) {
                     object@lowbnd <- value
                     return(object)
                 }
)


# uppbnd
setMethod("uppbnd", signature(object = "reactId_Exch"),
          function(object) {
              return(object@uppbnd)
          }
)

setReplaceMethod("uppbnd", signature = (object = "reactId_Exch"),
                 function(object, value) {
                     object@uppbnd <- value
                     return(object)
                 }
)


# uptake
setMethod("uptake", signature(object = "reactId_Exch"),
          function(object) {
              return(object@uptake)
          }
)

setReplaceMethod("uptake", signature = (object = "reactId_Exch"),
                 function(object, value) {
                     object@uptake <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("show", signature(object = "reactId_Exch"),
    function(object) {
        wcn <- trunc(log10(length(object))) + 3
        cn  <- paste("[", 1:length(object), "]", sep = "")
        lb <- ifelse(lowbnd(object) <= SYBIL_SETTINGS("MAXIMUM")*-1, -Inf, lowbnd(object))
        ub <- ifelse(uppbnd(object) >= SYBIL_SETTINGS("MAXIMUM"), Inf, uppbnd(object))
        cat(sprintf(paste("%", wcn, "s  ", "%-6s%-22s%-23s%-7s%6s %6s\n", sep = ""),
            "# ", "pos.", "reaction id", "metabolite id", "uptake", "lb", "ub"))
        cat(sprintf(paste("%", wcn, "s  ", "%-6s%-22s%-23s%-7s%6s %6s", sep = ""),
                    cn, react_pos(object), react_id(object),
                    met_id(object), uptake(object), lb, ub), sep = "\n")
        cat("\nnumber of exchange reactions: ", length(object), "\n", sep = "")
    }
)


# uptReact
setMethod("uptReact", signature(object = "reactId_Exch"),
          function(object) {
              return(object@react_id[object@uptake])
          }
)


# uptMet
setMethod("uptMet", signature(object = "reactId_Exch"),
          function(object) {
              return(object@met_id[object@uptake])
          }
)


# [
setMethod("[", signature(x = "reactId_Exch"),
    function(x, i, j, ..., drop = FALSE) {

        if ( (missing(i)) || (length(i) == 0) ) {
            return(x)
        }

        if (is(i, "character")) {
            ind <- match(i, react_id(x))
        }
#        else if (is(i, "reactId")) {
#            ind <- which(react_pos(x) %in% react_pos(i))
#        }
        else {
            ind <- i
        }

        if (max(ind, na.rm = TRUE) > length(x)) {
            stop("subscript out of bounds")
        }

        newRI <- new("reactId_Exch",
                     mod_id  = x@mod_id,
                     mod_key = x@mod_key,
                     rpnt    = x@react_pos[ind],
                     rid     = x@react_id[ind],
                     upt     = x@uptake[ind],
                     mpnt    = x@met_pos[ind],
                     mid     = x@met_id[ind],
                     lb      = x@lowbnd[ind],
                     ub      = x@uppbnd[ind])

        return(newRI)

    }
)
