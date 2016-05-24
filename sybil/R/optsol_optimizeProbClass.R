#  optsol_optimizeProbClass.R
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


# optsol_optimizeProbClass


#------------------------------------------------------------------------------#
#               definition of the class optsol_optimizeProb                    #
#------------------------------------------------------------------------------#

setClass("optsol_optimizeProb",
         representation(
              preProc  = "ppProc", # preprocessing lp result
              postProc = "ppProc"  # postprocessing lp result
        ),
        contains = "optsol"
)


#------------------------------------------------------------------------------#
#                              user constructor                                #
#------------------------------------------------------------------------------#

makeOptsolMO <- function(mod, sol) {

    # mod is of class modelorg
    # sol is the return value of optimizer()
    
    stopifnot(is(mod, "modelorg"), is(sol, "list"))

    opt <- new("optsol_optimizeProb",
        mod_id       = mod_id(mod),
        mod_key      = mod_key(mod),
        solver       = sol[["solver"]],
        method       = sol[["method"]],
        algorithm    = sol[["algorithm"]],
        num_of_prob  = as.integer(length(sol[["obj"]])),
        lp_num_cols  = as.integer(sol[["lp_num_cols"]]),
        lp_num_rows  = as.integer(sol[["lp_num_rows"]]),
        lp_obj       = as.numeric(sol[["obj"]]),
        lp_ok        = as.integer(sol[["ok"]]),
        lp_stat      = as.integer(sol[["stat"]]),
        lp_dir       = sol[["lp_dir"]],
        obj_coef     = obj_coef(mod),
        obj_func     = printObjFunc(mod),
        fldind       = as.integer(sol[["fldind"]]),
        fluxdist     = sol[["fluxdist"]],
        alg_par      = sol[["alg_par"]])


        if (!is.null(sol$prAna)) {
            preProc(opt) <- sol[["prAna"]]
        }

        if (!is.null(sol$poAna)) {
            postProc(opt) <- sol[["poAna"]]
        }

    return(opt)
}


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# preProc
setMethod("preProc", signature(object = "optsol_optimizeProb"),
          function(object) {
              return(object@preProc)
          }
)

setReplaceMethod("preProc", signature = (object = "optsol_optimizeProb"),
                 function(object, value) {
                     object@preProc <- value
                     return(object)
                 }
)


# postProc
setMethod("postProc", signature(object = "optsol_optimizeProb"),
          function(object) {
              return(object@postProc)
          }
)

setReplaceMethod("postProc", signature = (object = "optsol_optimizeProb"),
                 function(object, value) {
                     object@postProc <- value
                     return(object)
                 }
)
