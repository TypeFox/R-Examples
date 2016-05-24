#------------------------------------------------------------------------------#
#                             R interface to GLPK                              #
#------------------------------------------------------------------------------#

#  glpk_ptrClass.R
#  R interface to GLPK.
#
#  Copyright (C) 2011-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of glpkAPI.
#
#  GlpkAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GlpkAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with glpkAPI  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                          definition of class glpkPtr                         #
#------------------------------------------------------------------------------#


# representation of class glpkPtr
setClass(Class = "glpkPtr",
         representation(
              glpkPtrType = "character",
              glpkPointer = "externalptr"
#              ,
#              glpkParmSim = "externalptr",
#              glpkParmInt = "externalptr",
#              glpkParmMip = "externalptr"
         )
         #, contains = "externalptr"
)


#------------------------------------------------------------------------------#

# contructor for class glpkPtr
setMethod(f = "initialize",
          signature = "glpkPtr",
          definition = function(.Object, p, w) {

              .Object@glpkPointer <- attr(p, which = w, exact = TRUE)
              .Object@glpkPtrType <- as.character(p)

#              .Object@glpkParm    <- "parameter"
#              .Object@glpkParmSim <- attr(p, which = "simplex", exact = TRUE)
#              .Object@glpkParmInt <- attr(p, which = "interior", exact = TRUE)
#              .Object@glpkParmMip <- attr(p, which = "mip", exact = TRUE)
              
              return(.Object)
          
          }
)


# contructor for pointers to glpk problem structures
glpk_Pointer <- function(pointer) {

    if (is(pointer, "glpk_ptr")) {
        pObj <- new("glpkPtr",
                    p = pointer,
                    w = as.character("glpk_ptr"))
    }
    else {
        pObj <- pointer
    }

    return(pObj)
}

# contructor for pointers to translator workspace
trwks_Pointer <- function(pointer) {

    if (is(pointer, "trwks_ptr")) {
        pObj <- new("glpkPtr",
                    p = pointer,
                    w = as.character("trwks_ptr"))
    }
    else {
        pObj <- pointer
    }

    return(pObj)
}


#------------------------------------------------------------------------------#

# glpkPtrType
setMethod("glpkPtrType", signature(object = "glpkPtr"),
          function(object) {
              return(object@glpkPtrType)
          }
)

setReplaceMethod("glpkPtrType", signature = (object = "glpkPtr"),
                 function(object, value) {
                     object@glpkPtrType <- value
                     return(object)
                 }
)


# glpkPointer
setMethod("glpkPointer", signature(object = "glpkPtr"),
          function(object) {
              return(object@glpkPointer)
          }
)


#------------------------------------------------------------------------------#

setMethod("isNULLpointerGLPK", signature(object = "glpkPtr"),
    function(object) {
        return(.Call("isNULLptr", PACKAGE = "glpkAPI", glpkPointer(object)))
    }
)

setMethod("isGLPKpointer", signature(object = "glpkPtr"),
    function(object) {
        return(.Call("isGLPKptr", PACKAGE = "glpkAPI", glpkPointer(object)))
    }
)

setMethod("isTRWKSpointer", signature(object = "glpkPtr"),
    function(object) {
        return(.Call("isTRWKSptr", PACKAGE = "glpkAPI", glpkPointer(object)))
    }
)


#------------------------------------------------------------------------------#

setMethod("show", signature(object = "glpkPtr"),
    function(object) {
    
        nc <- NA
        
        if (isNULLpointerGLPK(object)) {
            ptrtype <- "NULL"
        }
        else {
            if (isGLPKpointer(object)) {
                ptrtype <- "GLPK problem object"
                nc <- getNumColsGLPK(object)
            }
            else if (isTRWKSpointer(object)) {
                ptrtype <- "MathProg translator workspace"
            }
            else {
                ptrtype <- "unknown"
            }
        }

        cat("object of class ", dQuote("glpkPtr"),
            ": pointer to ", ptrtype, ".\n", sep = "")

        if (!is.na(nc)) {
            if ( (nc < 1) || (nc > 10) ) {
                cat(paste("Number of variables:  ",
                          getNumColsGLPK(object), "\n"))
                cat(paste("Number of constraints:",
                          getNumRowsGLPK(object), "\n"))
            }
            else {
                # make a more illustrative method here
                cat(paste("Number of variables:  ",
                          getNumColsGLPK(object), "\n"))
                cat(paste("Number of constraints:",
                          getNumRowsGLPK(object), "\n"))
            }
        }
        
        cat(paste("Slot ",
                  dQuote("glpkPtrType"), ": ",
                  glpkPtrType(object), "\n", sep = ""))
        cat(paste("Slot ", dQuote("glpkPointer"), ": ", sep = ""))
        print(slot(object, "glpkPointer"), sep = "")
    }
)
