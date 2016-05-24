#------------------------------------------------------------------------------#
#                     R Interface to C API of IBM ILOG CPLEX                   #
#------------------------------------------------------------------------------#

#  cpxerrClass.R
#  R Interface to C API of IBM ILOG CPLEX Version 12.1 to 12.6.
#
#  Copyright (C) 2011-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of cplexAPI.
#
#  CplexAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  CplexAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with cplexAPI.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                       definition of class cplexError                         #
#------------------------------------------------------------------------------#


# representation of class cplexError
setClass(Class = "cplexError",
         representation(
              errnum = "integer"
         )
)


#------------------------------------------------------------------------------#

# contructor for class cplexError
cplexError <- function(err) {

    if(is(err, "cpxerr")) {
        cErr <- new("cplexError", errnum = as.integer(err))
    }
    else {
        cErr <- err
    }
    
    return(cErr)
}


#------------------------------------------------------------------------------#

# errnum
setMethod("errnum", signature(object = "cplexError"),
          function(object) {
              return(object@errnum)
          }
)

setReplaceMethod("errnum", signature = (object = "cplexError"),
                 function(object, value) {
                     object@errnum <- value
                     return(object)
                 }
)


# err
setMethod("err", signature(object = "cplexError"),
          function(object) {
              msg <- getErrorStrCPLEX(object@errnum)
              return(msg)
          }
)


# errmsg
setMethod("errmsg", signature(object = "cplexError"),
          function(object) {
              msg <- getErrorStrCPLEX(object@errnum)
              cat(msg)
          }
)


#------------------------------------------------------------------------------#

setMethod("show", signature(object = "cplexError"),
    function(object) {
    
        cat("object of class ", dQuote("cplexError"), ".\n", sep = "")
        cat(paste("Slot ", dQuote("errnum"), ": ",
                  errnum(object), "\n", sep = ""))
        cat(paste("Error string:  ", err(object), sep = ""))
    }
)
