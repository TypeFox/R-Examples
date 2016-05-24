#  optsol_fluxdelClass.R
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


# optsol_fluxdelClass


#------------------------------------------------------------------------------#
#                  definition of the class optsol_fluxdel                      #
#------------------------------------------------------------------------------#

setClass("optsol_fluxdel",
           representation(
               chlb      = "numeric",      # lower bound of changed fluxes/genes
               chub      = "numeric",      # upper bound of changed fluxes/genes
               dels      = "matrix"        # id's of deleted fluxes
           ),
           contains = "optsol_optimizeProb"
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# chlb
setMethod("chlb", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@chlb)
          }
)

setReplaceMethod("chlb", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@chlb <- value
                     return(object)
                 }
)


# chub
setMethod("chub", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@chub)
          }
)

setReplaceMethod("chub", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@chub <- value
                     return(object)
                 }
)


# dels
setMethod("dels", signature(object = "optsol_fluxdel"),
          function(object) {
              return(object@dels)
          }
)

setReplaceMethod("dels", signature = (object = "optsol_fluxdel"),
                 function(object, value) {
                     object@dels <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

# lethal
setMethod("lethal", signature(object = "optsol_fluxdel"),
          function(object, wt, tol) {
              
              stopifnot(is(wt, "numeric"), length(wt) == 1)

              if (missing(tol)) {
                  tol <- SYBIL_SETTINGS("TOLERANCE")
              }
              
              letid <- which(abs(mod_obj(object)/wt) < tol)
              let   <- logical(num_of_prob(object))
              let[letid] <- TRUE
              
              return(let)
          }
)


#setMethod("ind2id", signature = (object = "optsol_fluxdel"),
#                 function(object, slotN) {
#                     out <- NULL
#                     switch (slotN,
#                     
#                         "dels" = {
#                             out <- apply(dels(object), 2,
#                                          function(x) allGenes(object)[x]
#                                    )
#                         },
#                         
#                         {
#                             warning(paste("'", slotN, "' is not a valid slot!",
#                                           sep = ""
#                                    )
#                             )
#                         }
#                     
#                     )
#                 
#                     return(out)
#                 }
#)


setMethod("deleted", signature = (object = "optsol_fluxdel"),
                 function(object, i) {
                     value <- dels(object)[i, ]
                     return(value)
                 }
)


setMethod("[", "optsol_fluxdel", function(x, i, j, ..., drop = FALSE) {

        if ((missing(i)) || (length(i) == 0)) {
            return(x)
        }

        if (max(i) > length(x)) {
            stop("subscript out of bounds")
        }

        slots <- slotNames(x)
        
        isO <- is(x)[1]
        
        newSol <- new(isO,
            mod_id      = x@mod_id,
            mod_key     = x@mod_key,
            solver      = x@solver,
            method      = x@method,
            algorithm   = x@algorithm,
            num_of_prob = length(i),
            lp_num_cols = x@lp_num_cols,
            lp_num_rows = x@lp_num_rows,
            lp_obj      = x@lp_obj[i],
            lp_ok       = x@lp_ok[i],
            lp_stat     = x@lp_stat[i],
            lp_dir      = x@lp_dir,
            obj_coef    = x@obj_coef,
            obj_func    = x@obj_func,
            fldind      = x@fldind,
            chlb        = x@chlb[i],
            chub        = x@chub[i],
            dels        = x@dels[i, , drop = FALSE]
        )

        if (nfluxes(x) > 1) {
            NC_fl <- TRUE
        }
        else {
            NC_fl <- FALSE
        }

        if (isTRUE(NC_fl)) {
            newSol@fluxdist <- fluxDistribution(x@fluxes[ , i, drop = FALSE])
        }

        if ("fluxdels" %in% slots) {
            newSolfluxdels <- x@fluxdels[i]
        }

        if ("hasEffect" %in% slots) {
            newSol@hasEffect <- x@hasEffect[i]
        }

        return(newSol)


#         if ("num_of_prob" %in% slots) {
#             NC_num_of_prob <- length(i)
#             newClass <- paste(newClass, ", nprob = ", length(i)-1, sep = "")
#         }
#         else {
#             NC_num_prob <- NA
#         }
# 
#         if ("lp_dir" %in% slots) {
#             NC_lp_dir <- lp_dir(x)
#             newClass <- paste(newClass, ", lpdir = \"", lp_dir(x), "\"", sep = "")
#         }
#         else {
#             NC_lp_dir <- NA
#         }
# 
#         if ("lp_num_cols" %in% slots) {
#             NC_lp_num_cols <- lp_num_cols(x)
#             newClass <- paste(newClass, ", ncols = ", lp_num_cols(x), sep = "")
#         }
#         else {
#             NC_lp_num_cols <- NA
#         }
# 
#         if ("lp_num_rows" %in% slots) {
#             NC_lp_num_rows <- lp_num_rows(x)
#             newClass <- paste(newClass, ", nrows = ", lp_num_rows(x), sep = "")
#         }
#         else {
#             NC_lp_num_rows <- NA
#         }
# 
#         if ("obj_function" %in% slots) {
#             NC_obj_function <- obj_function(x)
#             newClass <- paste(newClass, ", objf = \"", obj_function(x), "\"", sep = "")
#         }
#         else {
#             NC_obj_function <- NA
#         }
# 
#         if ("fluxdist" %in% slots) {
#             if (is.na(fluxes(x))) {
#                 NC_fluxdist <- FALSE
#                 newClass <- paste(newClass, ", fld = ", FALSE, sep = "")
#             }
#             else {
#                 NC_fluxdist <- TRUE
#                 newClass <- paste(newClass, ", fld = ", TRUE, sep = "")
#             }
#         }
#         else {
#             NC_fluxdist <- NA
#         }
# 
#         if ("delmat" %in% slots) {
#             NC_delmat <- delmat(x)
#             dimdel <- dim(NC_delmat)
#             newClass <- paste(newClass, ", delrows = ", dimdel[1], ", delcols = ", dimdel[2], sep = "")
#         }
#         else {
#             NC_delmat <- NA
#         }
    }
)


