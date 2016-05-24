#  modelorgClass.R
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


# modelorgClass


#------------------------------------------------------------------------------#
#                      definition of the class modelorg                        #
#------------------------------------------------------------------------------#

# The class modelorg (the data part) is inspired by the data structure used
# in the COBRA Toolbox for the same purpose.

setClass("modelorg",
    representation(
         mod_desc     = "character",
         mod_name     = "character",   # model name
         mod_id       = "character",   # model id
         mod_key      = "character",   # model key (unique character string)
         mod_compart  = "character",   # vector compartments
         met_num      = "integer",     # number of metabolites
         met_id       = "character",   # vector metabolite id's
         met_name     = "character",   # vector metabolite names
         met_comp     = "integer",     # vector the metabolites compartment
         met_single   = "logical",     # metabolites appearing only once in S
         met_de       = "logical",     # dead end metabolites
         react_num    = "integer",     # number of reactions
         react_rev    = "logical",     # vector reversibilities
         react_id     = "character",   # vector reaction id's
         react_name   = "character",   # vector reaction names
         react_single = "logical",     # reactions using metabolites appearing only once in S
         react_de     = "logical",     # reactions using dead end metabolites
         S            = "Matrix",      # matrix S
         lowbnd       = "numeric",     # vector reactions lower bounds
         uppbnd       = "numeric",     # vector reactions upper bounds
         obj_coef     = "numeric",     # vector objective coefficients
         gprRules     = "character",
         genes        = "list",
         gpr          = "character",
         allGenes     = "character",
         rxnGeneMat   = "Matrix",
         subSys       = "Matrix"

    ),
    validity = .validmodelorg
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

modelorg <- function(id, name, subSys = NULL, compartment = NULL) {
    if (missing(id) || missing(name)) {
        stop("Creating an object of class model needs name and id!")
    }
    idP   <- as.character(id)
    nameP <- as.character(name)
    if (is.null(subSys)) {
        subSysP <- NULL
    }
    else {
        subSysP <- as.character(subSys)
    }
    if (is.null(compartment)) {
        compartmentP <- NULL
    }
    else {
        compartmentP <- as.character(compartment)
    }
    obj <- new("modelorg", id = idP, name = nameP,
               subSys = subSysP, compartment = compartmentP)
    return(obj)
}


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

setMethod(f = "initialize",
          signature = "modelorg",
          definition = function(.Object, id, name,
                                subSys = NULL, compartment = NULL) {

              if ( (!missing(id)) || (!missing(name)) ) {
                  .Object@mod_id     <- as.character(id)
                  .Object@mod_name   <- as.character(name)
                  .Object@mod_key    <- as.character(.generateModKey())
                  .Object@react_num  <- as.integer(0)
                  .Object@met_num    <- as.integer(0)
                  .Object@S          <- Matrix::Matrix(0, 0, 0)
                  .Object@rxnGeneMat <- Matrix::Matrix(0, 0, 0)
                  .Object@subSys     <- Matrix::Matrix(0, 0, length(subSys))
                  if (!is.null(subSys)) {
                      colnames(.Object@subSys) <- as.character(subSys)
                  }
                  if (!is.null(compartment)) {
                      .Object@mod_compart <- as.character(compartment)
                  }
              }

              return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# model id
setMethod("mod_id", signature(object = "modelorg"),
          function(object) {
              return(object@mod_id)
          }
)

setReplaceMethod("mod_id", signature = (object = "modelorg"),
                 function(object, value) {
                     object@mod_id <- value
                     return(object)
                 }
)


# model key
setMethod("mod_key", signature(object = "modelorg"),
          function(object) {
              return(object@mod_key)
          }
)

setReplaceMethod("mod_key", signature = (object = "modelorg"),
                 function(object, value) {
                     object@mod_key <- value
                     return(object)
                 }
)


# model name
setMethod("mod_name", signature(object = "modelorg"),
          function(object) {
              return(object@mod_name)
          }
)

setReplaceMethod("mod_name", signature = (object = "modelorg"),
                 function(object, value) {
                     object@mod_name <- value
                     return(object)
                 }
)


# model description
setMethod("mod_desc", signature(object = "modelorg"),
          function(object) {
              return(object@mod_desc)
          }
)

setReplaceMethod("mod_desc", signature = (object = "modelorg"),
                 function(object, value) {
                     object@mod_desc <- value
                     return(object)
                 }
)


# model compartments
setMethod("mod_compart", signature(object = "modelorg"),
          function(object) {
              return(object@mod_compart)
          }
)

setReplaceMethod("mod_compart", signature(object = "modelorg"),
          function(object, value) {
              object@mod_compart <- value
              return(object)
          }
)


# number of metabolites
setMethod("met_num", signature(object = "modelorg"),
          function(object) {
              return(object@met_num)
          }
)

setReplaceMethod("met_num", signature(object = "modelorg"),
          function(object, value) {
              object@met_num <- value
              return(object)
          }
)


# metabolite id's
setMethod("met_id", signature(object = "modelorg"),
          function(object) {
              return(object@met_id)
          }
)

setReplaceMethod("met_id", signature(object = "modelorg"),
          function(object, value) {
              object@met_id <- value
              return(object)
          }
)


# metabolite names
setMethod("met_name", signature(object = "modelorg"),
          function(object) {
              return(object@met_name)
          }
)

setReplaceMethod("met_name", signature(object = "modelorg"),
          function(object, value) {
              object@met_name <- value
              return(object)
          }
)


# metabolites compartments
setMethod("met_comp", signature(object = "modelorg"),
          function(object) {
              return(object@met_comp)
          }
)

setReplaceMethod("met_comp", signature(object = "modelorg"),
          function(object, value) {
              object@met_comp <- value
              return(object)
          }
)


# singletons
setMethod("met_single", signature(object = "modelorg"),
          function(object) {
              return(object@met_single)
          }
)

setReplaceMethod("met_single", signature(object = "modelorg"),
          function(object, value) {
              object@met_single <- value
              return(object)
          }
)


# dead ends
setMethod("met_de", signature(object = "modelorg"),
          function(object) {
              return(object@met_de)
          }
)

setReplaceMethod("met_de", signature(object = "modelorg"),
          function(object, value) {
              object@met_de <- value
              return(object)
          }
)


# number of reactions
setMethod("react_num", signature(object = "modelorg"),
          function(object) {
              return(object@react_num)
          }
)

setReplaceMethod("react_num", signature(object = "modelorg"),
          function(object, value) {
              object@react_num <- value
              return(object)
          }
)


# reversibilities
setMethod("react_rev", signature(object = "modelorg"),
          function(object) {
              return(object@react_rev)
          }
)

setReplaceMethod("react_rev", signature(object = "modelorg"),
          function(object, value) {
              object@react_rev <- value
              return(object)
          }
)


# reaction id's
setMethod("react_id", signature(object = "modelorg"),
          function(object) {
              return(object@react_id)
          }
)

setReplaceMethod("react_id", signature(object = "modelorg"),
          function(object, value) {
              object@react_id <- value
              return(object)
          }
)


# reaction names
setMethod("react_name", signature(object = "modelorg"),
          function(object) {
              return(object@react_name)
          }
)

setReplaceMethod("react_name", signature(object = "modelorg"),
          function(object, value) {
              object@react_name <- value
              return(object)
          }
)


# singletons
setMethod("react_single", signature(object = "modelorg"),
          function(object) {
              return(object@react_single)
          }
)

setReplaceMethod("react_single", signature(object = "modelorg"),
          function(object, value) {
              object@react_single <- value
              return(object)
          }
)


# dead ends
setMethod("react_de", signature(object = "modelorg"),
          function(object) {
              return(object@react_de)
          }
)

setReplaceMethod("react_de", signature(object = "modelorg"),
          function(object, value) {
              object@react_de <- value
              return(object)
          }
)


# stoichiometric matrix
setMethod("S", signature(object = "modelorg"),
          function(object) {
              return(object@S)
          }
)

setReplaceMethod("S", signature(object = "modelorg"),
          function(object, value) {
              object@S <- value
              return(object)
          }
)

# number of non zero elements in S
setMethod("Snnz", signature(object = "modelorg"),
          function(object) {
              return(length(object@S@x))
          }
)

# slot dimension of stoichiometric matrix (SparseM)
setMethod("dim", signature(x = "modelorg"),
          function(x) {
              return(dim(x@S))
          }
)


# lower bounds
setMethod("lowbnd", signature(object = "modelorg"),
          function(object) {
              return(object@lowbnd)
          }
)

setReplaceMethod("lowbnd", signature(object = "modelorg"),
          function(object, value) {
              object@lowbnd <- value
              return(object)
          }
)


# upper bounds
setMethod("uppbnd", signature(object = "modelorg"),
          function(object) {
              return(object@uppbnd)
          }
)

setReplaceMethod("uppbnd", signature(object = "modelorg"),
          function(object, value) {
              object@uppbnd <- value
              return(object)
          }
)


# objective coefficient
setMethod("obj_coef", signature(object = "modelorg"),
          function(object) {
              return(object@obj_coef)
          }
)

setReplaceMethod("obj_coef", signature(object = "modelorg"),
          function(object, value) {
              object@obj_coef <- value
              return(object)
          }
)


# gprRules
setMethod("gprRules", signature(object = "modelorg"),
          function(object) {
              return(object@gprRules)
          }
)

setReplaceMethod("gprRules", signature(object = "modelorg"),
          function(object, value) {
              object@gprRules <- value
              return(object)
          }
)


# genes
setMethod("genes", signature(object = "modelorg"),
          function(object) {
              return(object@genes)
          }
)

setReplaceMethod("genes", signature(object = "modelorg"),
          function(object, value) {
              object@genes <- value
              return(object)
          }
)


# gpr associations
setMethod("gpr", signature(object = "modelorg"),
          function(object) {
              return(object@gpr)
          }
)

setReplaceMethod("gpr", signature(object = "modelorg"),
          function(object, value) {
              object@gpr <- value
              return(object)
          }
)


# list of all genes
setMethod("allGenes", signature(object = "modelorg"),
          function(object) {
              return(object@allGenes)
          }
)

setReplaceMethod("allGenes", signature(object = "modelorg"),
          function(object, value) {
              object@allGenes <- value
              return(object)
          }
)


# reaction to gene mapping
setMethod("rxnGeneMat", signature(object = "modelorg"),
          function(object) {
              return(object@rxnGeneMat)
          }
)

setReplaceMethod("rxnGeneMat", signature(object = "modelorg"),
          function(object, value) {
              object@rxnGeneMat <- value
              return(object)
          }
)


# reaction sub systems
setMethod("subSys", signature(object = "modelorg"),
          function(object) {
              return(object@subSys)
          }
)

setReplaceMethod("subSys", signature(object = "modelorg"),
          function(object, value) {
              object@subSys <- value
              return(object)
          }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("show", signature(object = "modelorg"),
    function(object) {
        cat("model name:            ", mod_name(object), "\n")
        cat("number of compartments ", length(mod_compart(object)), "\n")
        lapply(mod_compart(object), function(x)
        cat("                       ", x, "\n")
              )
        cat("number of reactions:   ", react_num(object), "\n")
        cat("number of metabolites: ", met_num(object), "\n")
        if (length(allGenes) > 0) {
            cat("number of unique genes:", length(allGenes(object)), "\n")
        }
        cat("objective function:    ", printObjFunc(object), "\n")
    }
)


#------------------------------------------------------------------------------#

setMethod("optimizeProb", signature(object = "modelorg"),
    function(object, 
             algorithm = SYBIL_SETTINGS("ALGORITHM"),
             gene = NULL,
             react = NULL,
             lb = NULL,
             ub = NULL,
             retOptSol = TRUE,
             obj_coef = NULL,
             lpdir = NULL,
             mtfobj = NULL,
             fldind = TRUE,
             prCmd = NA, poCmd = NA,
             prCil = NA, poCil = NA, ...) {


        if (!is.null(gene)) {
            if (!is.null(react)) {
                warning("ignoring argument 'react'")
            }
            react <- geneDel(object, gene, checkId = TRUE)
            lb <- rep(lb[1], length(react))
            ub <- rep(ub[1], length(react))
        }

        # check the argument react
        # if object is of class "modelorg", react is given by the user
        if (!is.null(react)) {
            check <- checkReactId(object, react = react)
            if (is(check, "reactId")) {
                react <- react_pos(check)
                if (length(lb) == 1) {
                    lb <- rep(lb[1], length(react))
                }
                if (length(ub) == 1) {
                    ub <- rep(ub[1], length(react))
                }
                if (length(obj_coef) == 1) {
                    obj_coef <- rep(obj_coef[1], length(react))
                }
            }
            else {
                stop("check argument react")
            }
        }


        # -------------------------------------------------------------- #
        # run optimization
        # -------------------------------------------------------------- #

        if (algorithm == "mtf") {
            if (is.null(mtfobj)) {
                lpmod <- sysBiolAlg(model = object, algorithm = "mtf",
                                    react = react, lb = lb, ub = ub, ...)
            }
            else {
                lpmod <- sysBiolAlg(model = object,
                                    algorithm = "mtf", wtobj = mtfobj, ...)
            }
        }
        else {
            lpmod <- sysBiolAlg(model = object, algorithm = algorithm, ...)
        }
        
#        MoreArgs[["object"]]          <- lpmod
#        MoreArgs[["react"]]           <- react
#        MoreArgs[["lb"]]              <- lb
#        MoreArgs[["ub"]]              <- ub
#        
#        sol <- do.call("optimizeProb", args = MoreArgs)
        sol <- optimizeProb(lpmod,
                            react = react,
                            lb = lb,
                            ub = ub,
                            obj_coef = obj_coef,
                            lpdir = lpdir,
                            fldind = fldind,
                            resetChanges = FALSE,
                            prCmd = prCmd, poCmd = poCmd,
                            prCil = prCil, poCil = poCil)


        # -------------------------------------------------------------- #
        # store solution
        # -------------------------------------------------------------- #

        if (isTRUE(retOptSol)) {

            # solution object
            optsol <- new("optsol_optimizeProb",
                          mod_id       = mod_id(object),
                          mod_key      = mod_key(object),
                          solver       = solver(problem(lpmod)),
                          method       = method(problem(lpmod)),
                          algorithm    = algorithm(lpmod),
                          num_of_prob  = 1L,
                          lp_dir       = factor(getObjDir(problem(lpmod))),
                          lp_num_rows  = nr(lpmod),
                          lp_num_cols  = nc(lpmod),
                          lp_ok        = as.integer(sol[["ok"]]),
                          lp_obj       = sol$obj,
                          lp_stat      = as.integer(sol[["stat"]]),
                          obj_coef     = obj_coef(object),
                          obj_func     = printObjFunc(object),
                          fldind       = fldind(lpmod),
                          fluxdist     = fluxDistribution(fluxes = sol[["fluxes"]],
                                                    nrow = length(sol[["fluxes"]]),
                                                    ncol = 1L),
                          alg_par      = alg_par(lpmod))
 
            if (is(sol$preP, "ppProc")) {
                preProc(optsol) <- sol[["preP"]]
            }
    
            if (is(sol$postP, "ppProc")) {
                postProc(optsol) <- sol[["postP"]]
            }
    
            check <- validObject(optsol, test = TRUE)
    
            if (check != TRUE) {
                warning(paste("Validity check failed:", check, sep = "\n    "),
                        call. = FALSE
                )
            }
    
            checkOptSol(optsol, onlywarn = TRUE)
        }
        else {
            optsol             <- sol
            optsol[["fldind"]] <- fldind(lpmod)
        }

        delProb(problem(lpmod))
        remove(lpmod)
        
        return(optsol)

    }
)


#------------------------------------------------------------------------------#

# print objective function
setMethod("printObjFunc", signature(object = "modelorg"),
          function(object) {
              cInd <- obj_coef(object) != 0
              
              # check if there is an objective function
              if (sum(cInd) == 0) {
                  of <- "no objective function"
              }
              else {
                  obj <- gsub("^([^-])", "+\\1",
                              obj_coef(object)[cInd], perl = TRUE)
                  of  <- paste(paste(obj, react_id(object)[cInd]),
                               collapse = " ")
              }
            
              return(of)  
          }
)


#------------------------------------------------------------------------------#

# print reactions
setMethod("printReaction", signature(object = "modelorg"),
    function(object, react, printOut = TRUE, ...) {
  	
        check <- checkReactId(object, react = react)
        if (is(check, "reactId")) {
            cind <- react_pos(check)
        }
        else {
            stop("check argument react")
        }

        mat <- S(object)[, cind, drop = FALSE]
        nnz <- apply(mat, 2, "!=", 0)
        reaction <- character(length(cind))

        for (j in seq(along = cind)) {
        
            met <- met_id(object)[nnz[, j]]
            nzv <- mat[, j][nnz[, j]]
            
            ed <- nzv < 0
            pd <- nzv > 0

            if (sum(ed) > 0) {
                educt   <- paste(paste("(", abs(nzv[ed]), ")", sep = ""),
                                 met[ed], collapse = " + ")
            }
            else {
                educt = ""
            }

            if (sum(pd) > 0) {
                product <- paste(paste("(", nzv[pd], ")", sep = ""),
                                 met[pd], collapse = " + ")
            }
            else {
                product = ""
            }
            
            arrow   <- ifelse(react_rev(object)[cind[j]], " <==> ", " --> ")
            
            reaction[j] <- paste(react_id(check)[j],
                                 paste(educt, product, sep = arrow), sep = "\t")
        }

        if (isTRUE(printOut)) {
           cat("abbreviation\tequation", reaction, sep = "\n", ...)
        }
        
        return(invisible(reaction))

    }
)


#------------------------------------------------------------------------------#

# print metabolites
setMethod("printMetabolite", signature(object = "modelorg"),
    function(object, met, FBAlp = FALSE, printOut = TRUE, ...) {
  	
        if ( (missing(met)) || (isTRUE(FBAlp)) ) {
            mid  <- met_id(object)
            rind <- c(1:met_num(object))
        }
        else {
            if (is(met, "numeric")) {
                mid  <- met_id(object)[met]
                rind <- met
            }
            else if (is(met, "character")) {
                rind <- match(met, met_id(object))
                if (any(is.na(rind))) {
                    stop("check argument met")
                }
                else {
                    mid <- met
                }
            }
            else {
                stop("argument met must be character or numeric")
            }
        }
        
        # format metabolite id's to be compatible with CPLEX LP format
        midbr <- .makeLPcompatible(name = mid, prefix = "r", sep = "")

        # format reaction id's to be compatible with CPLEX LP format
        reaction_id <- .makeLPcompatible(name = react_id(object),
                                                 prefix = "x", sep = "")

        mat <- S(object)[rind, , drop = FALSE]
        nnz <- apply(mat, 1, "!=", 0)
        metabolite <- character(length(rind))


        for (i in seq(along = rind)) {
        
            react <- reaction_id[nnz[, i]]
            nzv   <- mat[i, ][nnz[, i]]
            
            nm <- nzv != 0

            if (sum(nm) > 0) {
                nz <- gsub("^([^-])", "+\\1", nzv[nm], perl = TRUE)
                bal <- paste(nz, react[nm], collapse = " ")
            }
            else {
                bal = ""
            }

            metabolite[i] <- paste(paste(" ", midbr[i], ":", sep = ""),
                                   bal, sep = "\t")
        }

        if (isTRUE(FBAlp)) {
            obj  <- c("Maximize",
                      paste(" obj:", printObjFunc(object)),
                      "",
                      "Subject To")
            bnds <- c("",
                      "Bounds",
                      paste(paste("", lowbnd(object)),
                            reaction_id, uppbnd(object), sep = " <= "),
                      "",
                      "End")
            
            tof  <- c(obj, paste(metabolite, "", sep = " = 0"), bnds)
        }
        else {
            tof  <- metabolite
        }
        
        if (isTRUE(printOut)) {
            cat(tof, sep = "\n", ...)
        }

        return(invisible(tof))

    }
)


#------------------------------------------------------------------------------#

setMethod("shrinkMatrix", signature(X = "modelorg"),
    function(X, i = NULL, j = NULL, tol = SYBIL_SETTINGS("TOLERANCE")) {
  	
        stopifnot(xor(is.null(i), is.null(j)))

        # look for reactions --> j
        if (is.null(i)) {
            # translate reaction id's to indices
            cj <- checkReactId(X, react = j)
            if (!is(cj, "reactId")) {
                stop("check argument j")
            }
            else {
                cn <- react_pos(cj)
            }
            met <- NULL
        }


        # look for metabolites --> i
        if (is.null(j)) {
            # translate reaction id's to indices
            if (is(i, "character")) {
                met <- which(met_id(X) %in% i)
            }
            else if (is(i, "reactId_Exch")) {
                met <- met_pos(i)
            }
            else if (is(i, "numeric")) {
                met <- i
            }
            else {
                stop("check argument i")
            }

            matB <- abs(S(X)[met, , drop = FALSE]) > tol

            cn <- colSums(matB) > 0
        }
         
        # generate new matrix
        mat <- S(X)[, cn, drop = FALSE]

        # binary matrix
        matB <- abs(mat) > tol

        if (is.null(met)) {
            rn <- rowSums(matB) > 0
        }
        else {
            nzm <- which(rowSums(matB) > 0)
            ord <- nzm %in% met
            rn <- c(nzm[ord], nzm[!ord])
        }

        mat <- mat[rn, , drop = FALSE]

        colnames(mat) <- react_id(X)[cn]
        rownames(mat) <- met_id(X)[rn]

        return(mat)
              
    }
)


#------------------------------------------------------------------------------#

setMethod("changeUptake", signature(object = "modelorg"),
    function(object, off = NULL, on = NULL,
             rate = SYBIL_SETTINGS("MAXIMUM") * -1) {
  	
        ex <- findExchReact(object)

        if (! is(off, "logical")) {
            if (is.null(off)) {
                reactOFF <- react_pos(ex)[uptake(ex)]
            }
            else {
                # metabolite id
                if (is(off, "character")) {
                    met      <- which(met_id(ex) %in% off)
                    reactOFF <- react_pos(ex)[met]
                }
                # metabolite index in S
                else if (is(off, "numeric")) {
                    met      <- which(met_pos(ex) %in% off)
                    reactOFF <- react_pos(ex)[met]
                }
                # exchange reactions
                else if (is(off, "reactId_Exch")) {
                    reactOFF <- react_pos(off)
                }
                else {
                    stop("check argument off")
                }
            
                if (length(reactOFF) < 1) {
                    stop("can not find ",
                         ngettext(is(off, "reactId_Exch"),
                                  "exchange reaction id ",
                                  "an exchange reaction for metabolite id "),
                         paste(sQuote(off), collapse = ", "))
                }
            }

            lowbnd(object)[reactOFF] <- 0
        }

        if (!is.null(on)) {
            # metabolite id
            if (is(on, "character")) {
                met     <- which(met_id(ex) %in% on)
                reactON <- react_pos(ex)[met]
            }
            # metabolite index in S
            else if (is(on, "numeric")) {
                met     <- which(met_pos(ex) %in% on)
                reactON <- react_pos(ex)[met]
            }
            # exchange reactions
            else if (is(on, "reactId_Exch")) {
                reactON <- react_pos(on)
            }
            else {
                stop("check argument on")
            }

            if (length(reactON) < 1) {
                stop("can not find ",
                     ngettext(is(on, "reactId_Exch"),
                              "exchange reaction id ",
                              "an exchange reaction for metabolite id "),
                     paste(sQuote(on), collapse = ", "))
            }

            lowbnd(object)[reactON] <- rate
        }
        
        return(object)
              
    }
)


#------------------------------------------------------------------------------#

setMethod("deadEndMetabolites", signature(object = "modelorg"),
    function(object, retIds = TRUE) {

  	    demr <- .deadEndMetabolite(mat = S(object), lb = lowbnd(object))
        
        if (isTRUE(retIds)) {
            dem <- list(dem = met_id(object)[demr[["dem"]]],
                        der = react_id(object)[demr[["der"]]])
        }
        else {
            dem <- demr
        }
        
        return(dem)
              
    }
)

#------------------------------------------------------------------------------#

setMethod("singletonMetabolites", signature(object = "modelorg"),
    function(object, tol = SYBIL_SETTINGS("TOLERANCE"), retIds = TRUE) {
  	
        Sb <- abs(S(object)) > tol

        singleton <- .singletonMetabolite(mat = Sb)

        if (isTRUE(retIds)) {
            sg <- list(smet = met_id(object)[singleton[["smet"]]],
                       sreact = react_id(object)[singleton[["sreact"]]])
        }
        else {
            sg <- singleton
        }

        return(sg)
              
    }
)














