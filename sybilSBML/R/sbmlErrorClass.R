#------------------------------------------------------------------------------#
#                          Link to libSBML for sybil                           #
#------------------------------------------------------------------------------#

#  sbmlErrorClass.R
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
#                         definition of class sbmlError                        #
#------------------------------------------------------------------------------#


# representation of class sbmlError
setClass(Class = "sbmlError",
         representation(
              sbmlInfos    = "list",
              sbmlWarnings = "list",
              sbmlErrors   = "list",
              sbmlFatals   = "list",
              sbmlFileName = "character",
              sbmlDocKey   = "character"
         )
)


#------------------------------------------------------------------------------#


# contructor function for class sbmlError
sbmlError <- function(err, sbmlf) {

    stopifnot(is(sbmlf, "sbmlPtr"))

    if (is(err, "sbml_error")) {
        pObj <- new("sbmlError",
                    sbmlInfos    = err[["infos"]],
                    sbmlWarnings = err[["warnings"]],
                    sbmlErrors   = err[["errors"]],
                    sbmlFatals   = err[["fatals"]],
                    sbmlFileName = sbmlFileName(sbmlf),
                    sbmlDocKey   = sbmlDocKey(sbmlf))
    }
    else {
        pObj <- err
    }

    return(pObj)
}


#------------------------------------------------------------------------------#

# sbmlInfos
setMethod("sbmlInfos", signature(object = "sbmlError"),
          function(object) {
              return(object@sbmlInfos)
          }
)

# sbmlWarnings
setMethod("sbmlWarnings", signature(object = "sbmlError"),
          function(object) {
              return(object@sbmlWarnings)
          }
)

# sbmlErrors
setMethod("sbmlErrors", signature(object = "sbmlError"),
          function(object) {
              return(object@sbmlErrors)
          }
)

# sbmlFatals
setMethod("sbmlFatals", signature(object = "sbmlError"),
          function(object) {
              return(object@sbmlFatals)
          }
)

# sbmlDocKey
setMethod("sbmlDocKey", signature(object = "sbmlError"),
          function(object) {
              return(object@sbmlDocKey)
          }
)

# sbmlFileName
setMethod("sbmlFileName", signature(object = "sbmlError"),
          function(object) {
              return(object@sbmlFileName)
          }
)


# getNumErrors
setMethod("getNumErrors", signature(object = "sbmlError"),
          function(object) {

              num    <- integer(5)
              num[1] <- length(sbmlInfos(object))
              num[2] <- length(sbmlWarnings(object))
              num[3] <- length(sbmlErrors(object))
              num[4] <- length(sbmlFatals(object))
              num[5] <- sum(num[1:4])
              names(num) <- c("Infos", "Warnings", "Errors", "Fatals", "Total")

              #cmd <- paste("length(sbml", ws, "(object))", sep = "")
              #num <- eval(parse(text = cmd))
              #names(num) <- ws

              return(num)
          }
)

#------------------------------------------------------------------------------#

# show
setMethod("show", signature(object = "sbmlError"),
    function(object) {

        cat("validation of SBML file ", sbmlFileName(object), "\n\n", sep = "")
        
        .printErrors(sbmlInfos(object), "Infos")
        .printErrors(sbmlWarnings(object), "Warnings")
        .printErrors(sbmlErrors(object), "Errors")
        .printErrors(sbmlFatals(object), "Fatals")
        
    }
)


# length
setMethod("length", signature(x = "sbmlError"),
          function(x) {
              num <- getNumErrors(x)
              names(num) <- NULL
              return(num[length(num)])
          }
)


setMethod("printSlot", signature(object = "sbmlError", ws = "character"),
    function(object, ws) {

        cmd <- paste(".printErrors(sbml", ws, "(object), '", ws, "')", sep = "")
        eval(parse(text = cmd))
        
    }
)


.printErrors <- function(err, type) {
    if (length(err) > 0) {
        cat(type, " (", length(err), "):\n", sep = "")
        i <- 0
        for (e in err) {
            i <- i + 1
            cat(sub("s$", "", type), " number ", i, ":\n", sep = "")
            cat("Id: ", e[["id"]], "\n", sep = "")
            cat("line: ", e[["line"]], ", column: ", e[["column"]], "\n", sep = "")
            cat("message:\n")
            cat(e[["message"]], "\n")
        }
    }
}
