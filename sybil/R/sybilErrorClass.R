#  sybilErrorClass.R
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


# sybilErrorClass


#------------------------------------------------------------------------------#
#                     definition of the class sybilError                       #
#------------------------------------------------------------------------------#

setClass("sybilError",
         representation(
              emsg = "character",
              enum = "integer"
         )
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

sybilError <- function(errmsg = "", number = NA) {

    msg <- paste(errmsg, collapse = " ")

    obj <- new(Class = "sybilError",
               emsg = as.character(msg),
               enum = as.integer(number)
    )

    return(obj)
}


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "sybilError",
          definition = function(.Object, emsg, enum) {

    if (missing(emsg)) {
        stop("creation of instances of class sybilError needs a message")
    }
    else {}

    if (missing(enum)) {
        enum <- NA
    }
    else {}

    .Object@emsg = as.character(emsg)
    .Object@enum = as.integer(enum)

    validObject(.Object)
    return(.Object)
}
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# msg
setMethod(f = "emsg",
          signature = "sybilError",
          definition = function(object) {
              return(object@emsg)
          }
)

setReplaceMethod(f = "emsg",
                 signature = "sybilError",
                 definition = function(object, value) {
                     object@emsg <- value
                     return(object)
                 }
)


# num
setMethod(f = "enum",
          signature = "sybilError",
          definition = function(object) {
              return(object@enum)
          }
)

setReplaceMethod(f = "enum",
                 signature = "sybilError",
                 definition = function(object, value) {
                     object@enum <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("show", signature(object = "sybilError"),
    function(object) {
        cat("error no.:", enum(object), "\n")
        cat(emsg(object), "\n")
    }
)

