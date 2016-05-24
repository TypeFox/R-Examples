#  modelorg_irrevClass.R
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


# modelorg_irrevClass


#------------------------------------------------------------------------------#
#                  definition of the class modelorg_irrev                      #
#------------------------------------------------------------------------------#

setClass("modelorg_irrev",
         representation(
                irrev     = "logical",
                matchrev  = "integer",
                rev2irrev = "matrix",
                irrev2rev = "integer"
         ),
         contains = "modelorg"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

modelorg_irrev <- function(id, name) {
    if (missing(id) || missing(name)) {
        stop("Creating an object of class modelorg_irrev needs name and id!")
    }
    id   <- as.character(id)
    name <- as.character(name)
    obj <- new("modelorg_irrev", id = id, name = name)
    return(obj)
}


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "modelorg_irrev",
          definition = function(.Object, id, name) {

              if (!missing(id) || !missing(name)) {
                  .Object <- callNextMethod(.Object, id = id, name = name)
              }
              
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# irrev
setMethod("irrev", signature(object = "modelorg_irrev"),
          function(object) {
              return(object@irrev)
          }
)

setReplaceMethod("irrev", signature = (object = "modelorg_irrev"),
                 function(object, value) {
                     object@irrev <- value
                     return(object)
                 }
)


# matchrev
setMethod("matchrev", signature(object = "modelorg_irrev"),
          function(object) {
              return(object@matchrev)
          }
)

setReplaceMethod("matchrev", signature = (object = "modelorg_irrev"),
                 function(object, value) {
                     object@matchrev <- value
                     return(object)
                 }
)


# rev2irrev
setMethod("rev2irrev", signature(object = "modelorg_irrev"),
         function(object) {
             return(object@rev2irrev)
         }
)

setReplaceMethod("rev2irrev", signature = (object = "modelorg_irrev"),
                function(object, value) {
                    object@rev2irrev <- value
                    return(object)
                }
)


# irrev2rev
setMethod("irrev2rev", signature(object = "modelorg_irrev"),
         function(object) {
             return(object@irrev2rev)
         }
)

setReplaceMethod("irrev2rev", signature = (object = "modelorg_irrev"),
                function(object, value) {
                    object@irrev2rev<- value
                    return(object)
                }
)



