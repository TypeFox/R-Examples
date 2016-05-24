#------------------------------------------------------------------------------#
#                     R Interface to C API of IBM ILOG CPLEX                   #
#------------------------------------------------------------------------------#

#  cplex_ptrClass.R
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
#                         definition of class cplexPtr                         #
#------------------------------------------------------------------------------#


# representation of class cplexPtr
setClass(Class = "cplexPtr",
         representation(
              cplexPtrType = "character",
              cplexPointer = "externalptr"
         )
         #, contains = "externalptr"
)


#------------------------------------------------------------------------------#

# contructor for class cplexPtr
setMethod(f = "initialize",
          signature = "cplexPtr",
          definition = function(.Object, p, w) {

              fn <- attr(p, which = "CPLEXfn", exact = TRUE)

              .Object@cplexPointer <- attr(p, which = w, exact = TRUE)
              .Object@cplexPtrType <- as.character(p)

              attr(.Object@cplexPtrType, which = "CPLEXfn") <- fn
              
              return(.Object)
          
          }
)


# contructor for pointers to cplex envoronments
cplex_ProbPointer <- function(pointer) {

    if(is(pointer, "cplex_prob_ptr")) {
        pObj <- new("cplexPtr",
                    p = pointer,
                    w = as.character("cplex_prob_ptr"))
    }
    else {
        pObj <- pointer
    }
    
    return(pObj)
}

# contructor for pointers to cplex problem onjects
cplex_EnvPointer <- function(pointer) {

    if(is(pointer, "cplex_env_ptr")) {
        pObj <- new("cplexPtr",
                    p = pointer,
                    w = as.character("cplex_env_ptr"))
    }
    else {
        pObj <- pointer
    }
    
    return(pObj)
}

# contructor for pointers to cplex files
cplex_FilePointer <- function(pointer) {

    if(is(pointer, "cplex_file_ptr")) {
        
        pObj <- new("cplexPtr",
                    p = pointer,
                    w = as.character("cplex_file_ptr"))
    }
    else {
        pObj <- pointer
    }
    
    return(pObj)
}

# contructor for pointers to cplex channels
cplex_ChannelPointer <- function(pointer, chname) {

    chn <- ifelse(missing(chname), "cplex_chan_ptr", chname)

    if(is(pointer, "cplex_chan_ptr")) {
        pObj <- new("cplexPtr",
                    p = pointer,
                    w = as.character(chn))
    }
    else {
        pObj <- pointer
    }
    
    return(pObj)
}

# contructor for pointers to termination signals
cplex_TermPointer <- function(pointer) {

    if(is(pointer, "cplex_term_ptr")) {
        pObj <- new("cplexPtr",
                    p = pointer,
                    w = as.character("cplex_term_ptr"))
    }
    else {
        pObj <- pointer
    }
    
    return(pObj)
}


#------------------------------------------------------------------------------#

setMethod("isNULLpointerCPLEX", signature(object = "cplexPtr"),
    function(object) {
        return(.Call("isNULLptr",
                     PACKAGE = "cplexAPI", cplexPointer(object)))
    }
)

setMethod("isCPLEXprobPointer", signature(object = "cplexPtr"),
    function(object) {
        return(.Call("isCPLEXprobPtr",
                     PACKAGE = "cplexAPI", cplexPointer(object)))
    }
)

setMethod("isCPLEXenvPointer", signature(object = "cplexPtr"),
    function(object) {
        return(.Call("isCPLEXenvPtr",
                     PACKAGE = "cplexAPI", cplexPointer(object)))
    }
)

setMethod("isCPLEXfilePointer", signature(object = "cplexPtr"),
    function(object) {
        return(.Call("isCPLEXfilePtr",
                     PACKAGE = "cplexAPI", cplexPointer(object)))
    }
)

setMethod("isCPLEXchanPointer", signature(object = "cplexPtr"),
    function(object) {
        return(.Call("isCPLEXchanPtr",
                     PACKAGE = "cplexAPI", cplexPointer(object)))
    }
)

setMethod("isCPLEXtermPointer", signature(object = "cplexPtr"),
    function(object) {
        return(.Call("isCPLEXtermPtr", 
                     PACKAGE = "cplexAPI", cplexPointer(object)))
    }
)


#------------------------------------------------------------------------------#

# cplexPtrType
setMethod("cplexPtrType", signature(object = "cplexPtr"),
          function(object) {
              return(object@cplexPtrType)
          }
)

setReplaceMethod("cplexPtrType", signature = (object = "cplexPtr"),
                 function(object, value) {
                     object@cplexPtrType <- value
                     return(object)
                 }
)


# cplexPointer
setMethod("cplexPointer", signature(object = "cplexPtr"),
          function(object) {
              return(object@cplexPointer)
          }
)


#------------------------------------------------------------------------------#

setMethod("show", signature(object = "cplexPtr"),
    function(object) {
    
        fn <- NA
        
        if (isNULLpointerCPLEX(object)) {
            ptrtype <- "NULL"
        }
        else {
            if (isCPLEXprobPointer(object)) {
                ptrtype <- "CPLEX problem object"
            }
            else if (isCPLEXenvPointer(object)) {
                ptrtype <- "CPLEX environment"
            }
            else if (isCPLEXfilePointer(object)) {
                ptrtype <- "CPLEX file"
                fn <- attr(cplexPtrType(object),
                           which = "CPLEXfn", exact = TRUE)
            }
            else if (isCPLEXchanPointer(object)) {
                ptrtype <- "CPLEX channel"
            }
            else if (isCPLEXtermPointer(object)) {
                ptrtype <- "CPLEX termination signal"
            }
            else {
                ptrtype <- "unknown"
            }
        }

        cat("object of class ", dQuote("cplexPtr"),
            ": pointer to ", ptrtype, ".\n", sep = "")
        cat(paste("Slot ",
                  dQuote("cplexPtrType"), ": ",
                  cplexPtrType(object), "\n", sep = ""))
        cat(paste("Slot ", dQuote("cplexPointer"), ": ", sep = ""))
        print(slot(object, "cplexPointer"), sep = "")
        if (!is.na(fn)) {
            cat(paste("Filename:   ", dQuote(fn), "\n"))
        }
    }
)


#------------------------------------------------------------------------------#

setMethod("summary", signature(object = "cplexPtr"),
    function(object, ...) {

        if (isNULLpointerCPLEX(object)) {
            cat("NULL pointer\n")
        }
        else {
            if (isCPLEXprobPointer(object)) {
                cat("CPLEX problem object\n")
                cat(paste("Number of variables:  ",
                          getNumColsCPLEX(lp = object, ...), "\n"))
                cat(paste("Number of constraints:",
                          getNumRowsCPLEX(lp = object, ...), "\n"))
                cat("\nSolution\n")
                sol <- solutionCPLEX(lp = object, ...)
                if (!is(sol, "cplexError")) {
                    nc <- getNumColsCPLEX(lp = object, ...)
                    nr <- getNumRowsCPLEX(lp = object, ...)
                    cat(paste("solution status:",
                              getStatStrCPLEX(stat = sol$lpstat, ...), "\n"))
                    
                    if (sol$lpstat == CPX_STAT_OPTIMAL) {
                        if (nc > 10) {
                            sx  <- paste(paste(round(sol$x[1:10], digits = 2),
                                               collapse = "  ")
                                         , "...")
                            sdj <- paste(paste(round(sol$dj[1:10], digits = 2),
                                               collapse = "  ")
                                         , "...")
                        }
                        else {
                            sx  <- paste(round(sol$x, digits = 2),
                                         collapse = "  ")
                            sdj <- paste(round(sol$dj, digits = 2),
                                         collapse = "  ")
                        }

                        if (nr > 10) {
                            pi <- paste(paste(round(sol$pi[1:10], digits = 2),
                                              collapse = "  ")
                                        , "...")
                            sl <- paste(paste(round(sol$slack[1:10],digits = 2),
                                              collapse = "  ")
                                        , "...")
                        }
                        else {
                            pi <- paste(round(sol$pi, digits = 2),
                                        collapse = "  ")
                            sl <- paste(round(sol$slack, digits = 2),
                                        collapse = "  ")
                        }

                        cat(paste("objective value:       ", sol$objval, "\n"))
                        cat(paste("variable values:       ", sx, "\n"))
                        cat(paste("dual variable values:  ", pi, "\n"))
                        cat(paste("slack variable values: ", sl, "\n"))
                        cat(paste("variable reduced costs:", sdj, "\n"))
                    }

                    cat("\nSolution information\n")
                    soln <- solnInfoCPLEX(lp = object, ...)
                    cat(paste("method:         ", soln$method, "\n"))
                    cat(paste("solution type:  ", soln$type, "\n"))
                    cat(paste("primal feasible:", soln$primal_feasible, "\n"))
                    cat(paste("dual feasible:  ", soln$dual_feasible, "\n"))

                }
            }
            else if (isCPLEXenvPointer(object)) {
                cat("CPLEX environment\n")
            }
            else if (isCPLEXfilePointer(object)) {
                cat("CPLEX file\n")
            }
            else if (isCPLEXchanPointer(object)) {
                cat("CPLEX channel\n")
            }
            else if (isCPLEXtermPointer(object)) {
                cat("CPLEX termination signal\n")
            }
            else {
                cat("unknown pointer\n")
            }
        }
        return(invisible(NULL))
    }
)
