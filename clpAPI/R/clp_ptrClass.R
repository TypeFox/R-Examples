#------------------------------------------------------------------------------#
#                           R interface to COIN-OR Clp                         #
#------------------------------------------------------------------------------#

#  clp_ptrClass.R
#  R interface to COIN-OR Clp.
#
#  Copyright (C) 2011-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of clpAPI.
#
#  ClpAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ClpAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with clpAPI  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                          definition of class clpPtr                          #
#------------------------------------------------------------------------------#


# representation of class clpPtr
setClass(Class = "clpPtr",
         representation(
              clpPtrType = "character",
              clpPointer = "externalptr"
         )
         #, contains = "externalptr"
)


#------------------------------------------------------------------------------#

# contructor for class clpPtr
setMethod(f = "initialize",
          signature = "clpPtr",
          definition = function(.Object, p, w) {

              .Object@clpPointer <- attr(p, which = w, exact = TRUE)
              .Object@clpPtrType <- as.character(p)
              
              return(.Object)
          
          }
)


# contructor for pointers to clp problem structures
clp_Pointer <- function(pointer) {

    if (is(pointer, "clp_ptr")) {
        pObj <- new("clpPtr",
                    p = pointer,
                    w = as.character("clp_ptr"))
    }
    else {
        pObj <- pointer
    }

    return(pObj)
}


#------------------------------------------------------------------------------#

# clpPtrType
setMethod("clpPtrType", signature(object = "clpPtr"),
          function(object) {
              return(object@clpPtrType)
          }
)

setReplaceMethod("clpPtrType", signature = (object = "clpPtr"),
                 function(object, value) {
                     object@clpPtrType <- value
                     return(object)
                 }
)


# clpPointer
setMethod("clpPointer", signature(object = "clpPtr"),
          function(object) {
              return(object@clpPointer)
          }
)


#------------------------------------------------------------------------------#

setMethod("isNULLpointerCLP", signature(object = "clpPtr"),
    function(object) {
        return(.Call("isNULLptr", PACKAGE = "clpAPI", clpPointer(object)))
    }
)

setMethod("isCLPpointer", signature(object = "clpPtr"),
    function(object) {
        return(.Call("isCLPptr", PACKAGE = "clpAPI", clpPointer(object)))
    }
)


#------------------------------------------------------------------------------#

setMethod("show", signature(object = "clpPtr"),
    function(object) {
    
        nc <- NA
        
        if (isNULLpointerCLP(object)) {
            ptrtype <- "NULL"
        }
        else {
            if (isCLPpointer(object)) {
                ptrtype <- "COIN-OR Clp problem object"
                nc <- getNumColsCLP(object)
            }
            else {
                ptrtype <- "unknown"
            }
        }

        cat("object of class ", dQuote("clpPtr"),
            ": pointer to ", ptrtype, ".\n", sep = "")

        if (!is.na(nc)) {
            if ( (nc < 1) || (nc > 10) ) {
                cat(paste("Number of variables:  ",
                          getNumColsCLP(object), "\n"))
                cat(paste("Number of constraints:",
                          getNumRowsCLP(object), "\n"))
            }
            else {
                # make a more illustrative method here
                cat(paste("Number of variables:  ",
                          getNumColsCLP(object), "\n"))
                cat(paste("Number of constraints:",
                          getNumRowsCLP(object), "\n"))
            }
        }
        
        cat(paste("Slot ",
                  dQuote("clpPtrType"), ": ",
                  clpPtrType(object), "\n", sep = ""))
        cat(paste("Slot ", dQuote("clpPointer"), ": ", sep = ""))
        print(slot(object, "clpPointer"), sep = "")
    }
)
