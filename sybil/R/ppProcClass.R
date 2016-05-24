#  ppProcClass.R
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


# ppProcClass


#------------------------------------------------------------------------------#
#                       definition of the class ppProc                         #
#------------------------------------------------------------------------------#

setClass("ppProc",
         representation(
              #cmd = "character",
              cmd = "list",
              pa  = "list",
              ind = "integer"
         )
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

ppProc <- function(cmd) {
    if (missing(cmd)) {
        stop("Creating an object of class ppProc needs a list of commands!")
    }

    if(is(cmd, "list")) {
        cmd <- as.list(cmd)
    }

    obj <- new("ppProc", cmd = cmd)
    return(obj)
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# cmd
setMethod("cmd", signature(object = "ppProc"),
          function(object) {
              return(object@cmd)
          }
)

setReplaceMethod("cmd", signature = (object = "ppProc"),
                 function(object, value) {
                     object@cmd <- value
                     return(object)
                 }
)


# pa
setMethod("pa", signature(object = "ppProc"),
          function(object) {
              return(object@pa)
          }
)

setReplaceMethod("pa", signature = (object = "ppProc"),
                 function(object, value) {
                     object@pa <- value
                     return(object)
                 }
)


# ind
setMethod("ind", signature(object = "ppProc"),
          function(object) {
              return(object@ind)
          }
)

setReplaceMethod("ind", signature = (object = "ppProc"),
                 function(object, value) {
                     object@ind <- value
                     return(object)
                 }
)
