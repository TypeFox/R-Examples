#------------------------------------------------------------------------------#
#                          Link to libSBML for sybil                           #
#------------------------------------------------------------------------------#

#  sbmlPtrClass.R
#  Link to libSBML for sybil.
#
#  Copyright (C) 2010-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybilSBML.
#
#  SybilSBML is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SybilSBML is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with SybilSBML.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                          definition of class sbmlPtr                         #
#------------------------------------------------------------------------------#


# representation of class sbmlPtr
setClass(Class = "sbmlPtr",
         representation(
              sbmlPtrType  = "character",
              sbmlPointer  = "externalptr",
              sbmlDocKey   = "character",
              sbmlFileName = "character"
         )
)


#------------------------------------------------------------------------------#

# contructor for class sbmlPtr
setMethod(f = "initialize",
          signature = "sbmlPtr",
          definition = function(.Object, p, w, key, fname) {

              .Object@sbmlPointer  <- attr(p, which = w, exact = TRUE)
              .Object@sbmlPtrType  <- as.character(p)
              .Object@sbmlDocKey   <- as.character(key)
              .Object@sbmlFileName <- as.character(fname)

              return(.Object)
          
          }
)


# contructor for pointers to sbml documents
sbmlDocPointer <- function(pointer) {

    if (is(pointer, "sbml_doc_ptr")) {
        pObj <- new("sbmlPtr",
                    p = pointer,
                    w = as.character("sbml_doc_ptr"),
                    key = as.character(sybil:::.generateModKey()),
                    fname = as.character(attr(pointer, which = "file_name", exact = TRUE)))
    }
    else {
        pObj <- pointer
    }

    return(pObj)
}


# contructor for pointers to sbml models
sbmlModPointer <- function(pointer, sbmlDoc) {

    if (is(pointer, "sbml_model_ptr")) {
        pObj <- new("sbmlPtr",
                    p = pointer,
                    w = as.character("sbml_model_ptr"),
                    key = sbmlDocKey(sbmlDoc),
                    fname = sbmlFileName(sbmlDoc))
    }
    else {
        pObj <- pointer
    }

    return(pObj)
}


#------------------------------------------------------------------------------#

# sbmlPtrType
setMethod("sbmlPtrType", signature(object = "sbmlPtr"),
          function(object) {
              return(object@sbmlPtrType)
          }
)

# sbmlPointer
setMethod("sbmlPointer", signature(object = "sbmlPtr"),
          function(object) {
              return(object@sbmlPointer)
          }
)

# sbmlDocKey
setMethod("sbmlDocKey", signature(object = "sbmlPtr"),
          function(object) {
              return(object@sbmlDocKey)
          }
)

# sbmlFileName
setMethod("sbmlFileName", signature(object = "sbmlPtr"),
          function(object) {
              return(object@sbmlFileName)
          }
)


#------------------------------------------------------------------------------#

setMethod("isNULLpointerSBML", signature(object = "sbmlPtr"),
    function(object) {
        return(.Call("isNULLptr", PACKAGE = "sybilSBML", sbmlPointer(object)))
    }
)

setMethod("isSBMLdocpointer", signature(object = "sbmlPtr"),
    function(object) {
        return(.Call("isSBMLdocptr", PACKAGE = "sybilSBML", sbmlPointer(object)))
    }
)

setMethod("isSBMLmodpointer", signature(object = "sbmlPtr"),
    function(object) {
        return(.Call("isSBMLmodptr", PACKAGE = "sybilSBML", sbmlPointer(object)))
    }
)


#------------------------------------------------------------------------------#

setMethod("show", signature(object = "sbmlPtr"),
    function(object) {
        
        if (isNULLpointerSBML(object)) {
            ptrtype <- "NULL"
        }
        else {
            if (isSBMLdocpointer(object)) {
                ptrtype <- "SBML document"
            }
            else if (isSBMLmodpointer(object)) {
                ptrtype <- "SBML model"
            }
            else {
                ptrtype <- "unknown"
            }
        }

        cat("object of class ", dQuote("sbmlPtr"),
            ": pointer to ", ptrtype, ".\n", sep = "")

        cat(paste("Slot ",
                  dQuote("sbmlPtrType"), ":  ",
                  sbmlPtrType(object), "\n", sep = ""))
        cat(paste("Slot ", dQuote("sbmlPointer"), ":  ", sep = ""))
        print(slot(object, "sbmlPointer"), sep = "")
        cat(paste("Slot ",
                  dQuote("sbmlDocKey"), ":   ",
                  sbmlDocKey(object), "\n", sep = ""))
        cat(paste("Slot ",
                  dQuote("sbmlFileName"), ": ",
                  sbmlFileName(object), "\n", sep = ""))
    }
)

