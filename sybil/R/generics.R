#  generics.R
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


#------------------------------------------------------------------------------#
#                                   generics                                   #
#------------------------------------------------------------------------------#

# from graphics
setGeneric("plot")

# sybil
setGeneric(name = "addCols",
           def  = function(lp, ncols, ...) { standardGeneric("addCols") }
)

setGeneric(name = "addColsToProb",
           def  = function(lp, ...) { standardGeneric("addColsToProb") }
)

setGeneric(name = "addRows",
           def  = function(lp, nrows, ...) { standardGeneric("addRows") }
)

setGeneric(name = "addRowsCols",
           def  = function(lp, nrows, ncols, ...) { standardGeneric("addRowsCols") }
)

setGeneric(name = "addRowsToProb",
           def  = function(lp, ...) { standardGeneric("addRowsToProb") }
)

setGeneric(name = "algorithm",
           def  = function(object) { standardGeneric("algorithm") }
)
setGeneric(name = "algorithm<-",
           def  = function(object, value) { standardGeneric("algorithm<-") }
)

setGeneric(name = "alg_par",
           def  = function(object) { standardGeneric("alg_par") }
)
setGeneric(name = "alg_par<-",
           def  = function(object, value) { standardGeneric("alg_par<-") }
)

setGeneric(name = "allGenes",
           def  = function(object) { standardGeneric("allGenes") }
)
setGeneric(name = "allGenes<-",
           def  = function(object, value) { standardGeneric("allGenes<-") }
)

setGeneric(name = "applyChanges",
           def  = function(object, del, obj, ld, react, lb, ub, obj_coef, fldind, lpdir) { standardGeneric("applyChanges") }
)

setGeneric(name = "backupProb",
           def  = function(lp) { standardGeneric("backupProb") }
)

setGeneric(name = "blocked",
           def  = function(object) { standardGeneric("blocked") }
)
setGeneric(name = "blocked<-",
           def  = function(object, value) { standardGeneric("blocked<-") }
)

setGeneric(name = "blReact",
           def  = function(object, tol = SYBIL_SETTINGS("TOLERANCE")) { standardGeneric("blReact") }
)

setGeneric(name = "changeColsBnds",
           def  = function(lp, ...) { standardGeneric("changeColsBnds") }
)

setGeneric(name = "changeColsBndsObjCoefs",
           def  = function(lp, ...) { standardGeneric("changeColsBndsObjCoefs") }
)

setGeneric(name = "changeMatrixRow",
           def  = function(lp, i, j, val) { standardGeneric("changeMatrixRow") }
)

setGeneric(name = "changeMaxObj",
           def  = function(object, ...) { standardGeneric("changeMaxObj") }
)

setGeneric(name = "changeObjCoefs",
           def  = function(lp, ...) { standardGeneric("changeObjCoefs") }
)

setGeneric(name = "changeRowsBnds",
           def  = function(lp, ...) { standardGeneric("changeRowsBnds") }
)

setGeneric(name = "changeUptake",
           def  = function(object, ...) { standardGeneric("changeUptake") }
)

setGeneric(name = "checkOptSol",
           def  = function(opt, ...) { standardGeneric("checkOptSol") }
)

setGeneric(name = "checkStat",
           def  = function(opt) { standardGeneric("checkStat") }
)

setGeneric(name = "chlb",
           def  = function(object) { standardGeneric("chlb") }
)
setGeneric(name = "chlb<-",
           def  = function(object, value) { standardGeneric("chlb<-") }
)

setGeneric(name = "chub",
           def  = function(object) { standardGeneric("chub") }
)
setGeneric(name = "chub<-",
           def  = function(object, value) { standardGeneric("chub<-") }
)

setGeneric(name = "cmd",
           def  = function(object) { standardGeneric("cmd") }
)
setGeneric(name = "cmd<-",
           def  = function(object, value) { standardGeneric("cmd<-") }
)

setGeneric(name = "ctrlfl",
           def  = function(object) { standardGeneric("ctrlfl") }
)
setGeneric(name = "ctrlfl<-",
           def  = function(object, value) { standardGeneric("ctrlfl<-") }
)

setGeneric(name = "ctrlr",
           def  = function(object) { standardGeneric("ctrlr") }
)
setGeneric(name = "ctrlr<-",
           def  = function(object, value) { standardGeneric("ctrlr<-") }
)

setGeneric(name = "deleted",
           def  = function(object, ...) { standardGeneric("deleted") }
)

setGeneric(name = "deadEndMetabolites",
           def  = function(object, ...) { standardGeneric("deadEndMetabolites") }
)
# setGeneric(name = "deleted<-",
#            def  = function(object, value) { standardGeneric("deleted<-") }
# )

setGeneric(name = "delProb",
           def  = function(lp, ...) { standardGeneric("delProb") }
)

setGeneric(name = "dels",
           def  = function(object) { standardGeneric("dels") }
)
setGeneric(name = "dels<-",
           def  = function(object, value) { standardGeneric("dels<-") }
)

setGeneric(name = "didFoot",
           def  = function(object) { standardGeneric("didFoot") }
)
setGeneric(name = "didFoot<-",
           def  = function(object, value) { standardGeneric("didFoot<-") }
)

setGeneric(name = "emsg",
           def  = function(object) { standardGeneric("emsg") }
)
setGeneric(name = "emsg<-",
           def  = function(object, value) { standardGeneric("emsg<-") }
)

setGeneric(name = "enum",
           def  = function(object) { standardGeneric("enum") }
)
setGeneric(name = "enum<-",
           def  = function(object, value) { standardGeneric("enum<-") }
)

setGeneric(name = "exit_code",
           def  = function(object) { standardGeneric("exit_code") }
)
setGeneric(name = "exit_code<-",
           def  = function(object, value) { standardGeneric("exit_code<-") }
)

setGeneric(name = "exit_meaning",
           def  = function(object) { standardGeneric("exit_meaning") }
)
setGeneric(name = "exit_meaning<-",
           def  = function(object, value) { standardGeneric("exit_meaning<-") }
)

setGeneric(name = "exit_num",
           def  = function(object) { standardGeneric("exit_num") }
)
setGeneric(name = "exit_num<-",
           def  = function(object, value) { standardGeneric("exit_num<-") }
)

setGeneric(name = "ex_met",
           def  = function(object) { standardGeneric("ex_met") }
)

setGeneric(name = "ex_val",
           def  = function(object) { standardGeneric("ex_val") }
)

setGeneric(name = "fenc",
           def  = function(object) { standardGeneric("fenc") }
)
setGeneric(name = "fenc<-",
           def  = function(object, value) { standardGeneric("fenc<-") }
)

setGeneric(name = "fh",
           def  = function(object) { standardGeneric("fh") }
)
setGeneric(name = "fh<-",
           def  = function(object, value) { standardGeneric("fh<-") }
)

setGeneric(name = "fldind",
           def  = function(object) { standardGeneric("fldind") }
)
setGeneric(name = "fldind<-",
           def  = function(object, value) { standardGeneric("fldind<-") }
)

setGeneric(name = "fluxdels",
           def  = function(object) { standardGeneric("fluxdels") }
)
setGeneric(name = "fluxdels<-",
           def  = function(object, value) { standardGeneric("fluxdels<-") }
)

setGeneric(name = "fluxdist",
           def  = function(object) { standardGeneric("fluxdist") }
)
setGeneric(name = "fluxdist<-",
           def  = function(object, value) { standardGeneric("fluxdist<-") }
)

setGeneric(name = "fluxes",
           def  = function(object) { standardGeneric("fluxes") }
)
setGeneric(name = "fluxes<-",
           def  = function(object, value) { standardGeneric("fluxes<-") }
)

setGeneric(name = "fname",
           def  = function(object) { standardGeneric("fname") }
)
setGeneric(name = "fname<-",
           def  = function(object, value) { standardGeneric("fname<-") }
)

setGeneric(name = "fpath",
           def  = function(object) { standardGeneric("fpath") }
)
setGeneric(name = "fpath<-",
           def  = function(object, value) { standardGeneric("fpath<-") }
)

setGeneric(name = "genes",
           def  = function(object) { standardGeneric("genes") }
)
setGeneric(name = "genes<-",
           def  = function(object, value) { standardGeneric("genes<-") }
)

setGeneric(name = "getColPrim",
           def  = function(lp, j) { standardGeneric("getColPrim") }
)

setGeneric(name = "getColsLowBnds",
           def  = function(lp, j) { standardGeneric("getColsLowBnds") }
)

setGeneric(name = "getColsNames",
           def  = function(lp, j) { standardGeneric("getColsNames") }
)

setGeneric(name = "getColsUppBnds",
           def  = function(lp, j) { standardGeneric("getColsUppBnds") }
)

setGeneric(name = "getFluxDist",
           def  = function(lp, ...) { standardGeneric("getFluxDist") }
)

setGeneric(name = "getNumCols",
           def  = function(lp) { standardGeneric("getNumCols") }
)

setGeneric(name = "getNumNnz",
           def  = function(lp) { standardGeneric("getNumNnz") }
)

setGeneric(name = "getNumRows",
           def  = function(lp) { standardGeneric("getNumRows") }
)

setGeneric(name = "getObjCoefs",
           def  = function(lp, j) { standardGeneric("getObjCoefs") }
)

setGeneric(name = "getObjDir",
           def  = function(lp) { standardGeneric("getObjDir") }
)

setGeneric(name = "getObjVal",
           def  = function(lp) { standardGeneric("getObjVal") }
)

setGeneric(name = "getRedCosts",
           def  = function(lp) { standardGeneric("getRedCosts") }
)

setGeneric(name = "getRowsLowBnds",
           def  = function(lp, i) { standardGeneric("getRowsLowBnds") }
)

setGeneric(name = "getRowsNames",
           def  = function(lp, i) { standardGeneric("getRowsNames") }
)

setGeneric(name = "getRowsUppBnds",
           def  = function(lp, i) { standardGeneric("getRowsUppBnds") }
)

setGeneric(name = "getSolStat",
           def  = function(lp) { standardGeneric("getSolStat") }
)

setGeneric(name = "getSolverParm",
           def  = function(lp) { standardGeneric("getSolverParm") }
)

setGeneric(name = "gpr",
           def  = function(object) { standardGeneric("gpr") }
)
setGeneric(name = "gpr<-",
           def  = function(object, value) { standardGeneric("gpr<-") }
)

setGeneric(name = "gprRules",
           def  = function(object) { standardGeneric("gprRules") }
)
setGeneric(name = "gprRules<-",
           def  = function(object, value) { standardGeneric("gprRules<-") }
)

setGeneric(name = "hasEffect",
           def  = function(object) { standardGeneric("hasEffect") }
)
setGeneric(name = "hasEffect<-",
           def  = function(object, value) { standardGeneric("hasEffect<-") }
)

setGeneric(name = "ind",
           def  = function(object) { standardGeneric("ind") }
)
setGeneric(name = "ind<-",
           def  = function(object, value) { standardGeneric("ind<-") }
)

#setGeneric(name = "ind2id",
#           def  = function(object, ...) { standardGeneric("ind2id") }
#)

setGeneric(name = "initProb",
           def  = function(lp, ...) { standardGeneric("initProb") }
)

setGeneric(name = "irrev",
           def  = function(object) { standardGeneric("irrev") }
)
setGeneric(name = "irrev<-",
           def  = function(object, value) { standardGeneric("irrev<-") }
)

setGeneric(name = "irrev2rev",
           def  = function(object) { standardGeneric("irrev2rev") }
)
setGeneric(name = "irrev2rev<-",
           def  = function(object, value) { standardGeneric("irrev2rev<-") }
)

setGeneric(name = "lethal",
           def  = function(object, wt, tol) { standardGeneric("lethal") }
)

setGeneric(name = "loadLPprob",
           def  = function(lp, ...) { standardGeneric("loadLPprob") }
)

setGeneric(name = "loadQobj",
           def  = function(lp, mat) { standardGeneric("loadQobj") }
)

setGeneric(name = "logCall",
           def  = function(object, nog) { standardGeneric("logCall") }
)
# setGeneric(name = "logCall",
#            def  = function(object, func, fargl, thdargs) { standardGeneric("logCall") }
# )
setGeneric(name = "logClose<-",
           def  = function(object, value) { standardGeneric("logClose<-") }
)
setGeneric(name = "logComment",
           def  = function(object, cmt, cmtChar) { standardGeneric("logComment") }
)
setGeneric(name = "logError",
           def  = function(object, msg, num) { standardGeneric("logError") }
)
setGeneric(name = "logFoot<-",
           def  = function(object, value) { standardGeneric("logFoot<-") }
)
setGeneric(name = "logFH",
           def  = function(object) { standardGeneric("logFH") }
)
setGeneric(name = "logHead",
           def  = function(object) { standardGeneric("logHead") }
)
setGeneric(name = "logMessage",
           def  = function(object, appendEllipsis, ...) { standardGeneric("logMessage") }
)
setGeneric(name = "logOptimization",
           def  = function(object, ...) { standardGeneric("logOptimization") }
)
setGeneric(name = "logOptimizationTH",
           def  = function(object) { standardGeneric("logOptimizationTH") }
)
setGeneric(name = "logStep<-",
           def  = function(object, value) { standardGeneric("logStep<-") }
)
setGeneric(name = "lstname",
           def  = function(object) { standardGeneric("lstname") }
)
setGeneric(name = "logWarning",
           def  = function(object, ...) { standardGeneric("logWarning") }
)

setGeneric(name = "loglevel",
           def  = function(object) { standardGeneric("loglevel") }
)
setGeneric(name = "loglevel<-",
           def  = function(object, value) { standardGeneric("loglevel<-") }
)

setGeneric(name = "lowbnd",
           def  = function(object) { standardGeneric("lowbnd") }
)
setGeneric(name = "lowbnd<-",
           def  = function(object, value) { standardGeneric("lowbnd<-") }
)

setGeneric(name = "lp_dir",
           def  = function(object) { standardGeneric("lp_dir") }
)
setGeneric(name = "lp_dir<-",
           def  = function(object, value) { standardGeneric("lp_dir<-") }
)

setGeneric(name = "lp_num_cols",
           def  = function(object) { standardGeneric("lp_num_cols") }
)
setGeneric(name = "lp_num_cols<-",
           def  = function(object, value) { standardGeneric("lp_num_cols<-") }
)

setGeneric(name = "lp_num_rows",
           def  = function(object) { standardGeneric("lp_num_rows") }
)
setGeneric(name = "lp_num_rows<-",
           def  = function(object, value) { standardGeneric("lp_num_rows<-") }
)

setGeneric(name = "lp_obj",
           def  = function(object) { standardGeneric("lp_obj") }
)
setGeneric(name = "lp_obj<-",
           def  = function(object, value) { standardGeneric("lp_obj<-") }
)

setGeneric(name = "lp_ok",
           def  = function(object) { standardGeneric("lp_ok") }
)
setGeneric(name = "lp_ok<-",
           def  = function(object, value) { standardGeneric("lp_ok<-") }
)

setGeneric(name = "lp_stat",
           def  = function(object) { standardGeneric("lp_stat") }
)
setGeneric(name = "lp_stat<-",
           def  = function(object, value) { standardGeneric("lp_stat<-") }
)

setGeneric(name = "matchrev",
           def  = function(object) { standardGeneric("matchrev") }
)
setGeneric(name = "matchrev<-",
           def  = function(object, value) { standardGeneric("matchrev<-") }
)

setGeneric(name = "maxSol",
           def  = function(object, ...) { standardGeneric("maxSol") }
)

setGeneric(name = "met_comp",
           def  = function(object) { standardGeneric("met_comp") }
)
setGeneric(name = "met_comp<-",
           def  = function(object, value) { standardGeneric("met_comp<-") }
)

setGeneric(name = "met_de",
           def  = function(object) { standardGeneric("met_de") }
)
setGeneric(name = "met_de<-",
           def  = function(object, value) { standardGeneric("met_de<-") }
)

setGeneric(name = "met_id",
           def  = function(object) { standardGeneric("met_id") }
)
setGeneric(name = "met_id<-",
           def  = function(object, value) { standardGeneric("met_id<-") }
)

setGeneric(name = "met_name",
           def  = function(object) { standardGeneric("met_name") }
)
setGeneric(name = "met_name<-",
           def  = function(object, value) { standardGeneric("met_name<-") }
)

setGeneric(name = "met_num",
           def  = function(object) { standardGeneric("met_num") }
)
setGeneric(name = "met_num<-",
           def  = function(object, value) { standardGeneric("met_num<-") }
)

setGeneric(name = "met_pos",
           def  = function(object) { standardGeneric("met_pos") }
)
setGeneric(name = "met_pos<-",
           def  = function(object, value) { standardGeneric("met_pos<-") }
)

setGeneric(name = "met_single",
           def  = function(object) { standardGeneric("met_single") }
)
setGeneric(name = "met_single<-",
           def  = function(object, value) { standardGeneric("met_single<-") }
)

setGeneric(name = "method",
           def  = function(object) { standardGeneric("method") }
)
setGeneric(name = "method<-",
           def  = function(object, value) { standardGeneric("method<-") }
)

setGeneric(name = "minSol",
           def  = function(object, ...) { standardGeneric("minSol") }
)

setGeneric(name = "mod_compart",
           def  = function(object) { standardGeneric("mod_compart") }
)
setGeneric(name = "mod_compart<-",
           def  = function(object, value) { standardGeneric("mod_compart<-") }
)

setGeneric(name = "mod_desc",
           def  = function(object) { standardGeneric("mod_desc") }
)
setGeneric(name = "mod_desc<-",
           def  = function(object, value) { standardGeneric("mod_desc<-") }
)

setGeneric(name = "mod_id",
           def  = function(object) { standardGeneric("mod_id") }
)
setGeneric(name = "mod_id<-",
           def  = function(object, value) { standardGeneric("mod_id<-") }
)

setGeneric(name = "mod_key",
           def  = function(object) { standardGeneric("mod_key") }
)
setGeneric(name = "mod_key<-",
           def  = function(object, value) { standardGeneric("mod_key<-") }
)

setGeneric(name = "mod_name",
           def  = function(object) { standardGeneric("mod_name") }
)
setGeneric(name = "mod_name<-",
           def  = function(object, value) { standardGeneric("mod_name<-") }
)

setGeneric(name = "mod_obj",
           def  = function(object) { standardGeneric("mod_obj") }
)
setGeneric(name = "mod_obj<-",
           def  = function(object, value) { standardGeneric("mod_obj<-") }
)

setGeneric(name = "nc",
           def  = function(object) { standardGeneric("nc") }
)
setGeneric(name = "nc<-",
           def  = function(object, value) { standardGeneric("nc<-") }
)

setGeneric(name = "nfluxes",
           def  = function(object) { standardGeneric("nfluxes") }
)

setGeneric(name = "nr",
           def  = function(object) { standardGeneric("nr") }
)
setGeneric(name = "nr<-",
           def  = function(object, value) { standardGeneric("nr<-") }
)

setGeneric(name = "nvar",
           def  = function(object) { standardGeneric("nvar") }
)

setGeneric(name = "nzeros",
           def  = function(object) { standardGeneric("nzeros") }
)

setGeneric(name = "num_of_fluxes",
           def  = function(object) { standardGeneric("num_of_fluxes") }
)

setGeneric(name = "num_of_prob",
           def  = function(object) { standardGeneric("num_of_prob") }
)
setGeneric(name = "num_of_prob<-",
           def  = function(object, value) { standardGeneric("num_of_prob<-") }
)

setGeneric(name = "obj_coef",
           def  = function(object) { standardGeneric("obj_coef") }
)
setGeneric(name = "obj_coef<-",
           def  = function(object, value) { standardGeneric("obj_coef<-") }
)

setGeneric(name = "obj_func",
           def  = function(object) { standardGeneric("obj_func") }
)
setGeneric(name = "obj_func<-",
           def  = function(object, value) { standardGeneric("obj_func<-") }
)

setGeneric(name = "optimizeProb",
           def  = function(object, ...) { standardGeneric("optimizeProb") }
)

setGeneric(name = "pa",
           def  = function(object) { standardGeneric("pa") }
)
setGeneric(name = "pa<-",
           def  = function(object, value) { standardGeneric("pa<-") }
)

setGeneric(name = "postProc",
           def  = function(object) { standardGeneric("postProc") }
)
setGeneric(name = "postProc<-",
           def  = function(object, value) { standardGeneric("postProc<-") }
)

setGeneric(name = "plotRangeVar",
           def  = function(object, ...) { standardGeneric("plotRangeVar") }
)

setGeneric(name = "preProc",
           def  = function(object) { standardGeneric("preProc") }
)
setGeneric(name = "preProc<-",
           def  = function(object, value) { standardGeneric("preProc<-") }
)

setGeneric(name = "printExchange",
           def  = function(object, ...) { standardGeneric("printExchange") }
)

setGeneric(name = "printMetabolite",
           def  = function(object, ...) { standardGeneric("printMetabolite") }
)

setGeneric(name = "printObjFunc",
           def  = function(object) { standardGeneric("printObjFunc") }
)

setGeneric(name = "printReaction",
           def  = function(object, ...) { standardGeneric("printReaction") }
)
setGeneric(name = "printReaction",
           def  = function(object, mod, ...) { standardGeneric("printReaction") }
)

setGeneric(name = "problem",
           def  = function(object) { standardGeneric("problem") }
)
#setGeneric(name = "problem<-",
#           def  = function(object, value) { standardGeneric("problem<-") }
#)

setGeneric(name = "probType",
           def  = function(object) { standardGeneric("probType") }
)

setGeneric(name = "rate",
           def  = function(object) { standardGeneric("rate") }
)

setGeneric(name = "react",
           def  = function(object) { standardGeneric("react") }
)
setGeneric(name = "react<-",
           def  = function(object, value) { standardGeneric("react<-") }
)

setGeneric(name = "react_de",
           def  = function(object) { standardGeneric("react_de") }
)
setGeneric(name = "react_de<-",
           def  = function(object, value) { standardGeneric("react_de<-") }
)

setGeneric(name = "react_id",
           def  = function(object) { standardGeneric("react_id") }
)
setGeneric(name = "react_id<-",
           def  = function(object, value) { standardGeneric("react_id<-") }
)

setGeneric(name = "react_name",
           def  = function(object) { standardGeneric("react_name") }
)
setGeneric(name = "react_name<-",
           def  = function(object, value) { standardGeneric("react_name<-") }
)

setGeneric(name = "react_num",
           def  = function(object) { standardGeneric("react_num") }
)
setGeneric(name = "react_num<-",
           def  = function(object, value) { standardGeneric("react_num<-") }
)

setGeneric(name = "react_pos",
           def  = function(object) { standardGeneric("react_pos") }
)
setGeneric(name = "react_pos<-",
           def  = function(object, value) { standardGeneric("react_pos<-") }
)

setGeneric(name = "react_rev",
           def  = function(object) { standardGeneric("react_rev") }
)
setGeneric(name = "react_rev<-",
           def  = function(object, value) { standardGeneric("react_rev<-") }
)

setGeneric(name = "react_single",
           def  = function(object) { standardGeneric("react_single") }
)
setGeneric(name = "react_single<-",
           def  = function(object, value) { standardGeneric("react_single<-") }
)

setGeneric(name = "readProb",
           def  = function(lp, fname, ff = "lp", ...) { standardGeneric("readProb") }
)

setGeneric(name = "resetChanges",
           def  = function(object, old_val) { standardGeneric("resetChanges") }
)

setGeneric(name = "rev2irrev",
           def  = function(object) { standardGeneric("rev2irrev") }
)
setGeneric(name = "rev2irrev<-",
           def  = function(object, value) { standardGeneric("rev2irrev<-") }
)

setGeneric(name = "rxnGeneMat",
           def  = function(object) { standardGeneric("rxnGeneMat") }
)
setGeneric(name = "rxnGeneMat<-",
           def  = function(object, value) { standardGeneric("rxnGeneMat<-") }
)

setGeneric(name = "S",
           def  = function(object) { standardGeneric("S") }
)
setGeneric(name = "S<-",
           def  = function(object, value) { standardGeneric("S<-") }
)

setGeneric(name = "scaleProb",
           def  = function(lp, ...) { standardGeneric("scaleProb") }
)

setGeneric(name = "shrinkMatrix",
           def  = function(X, ...) { standardGeneric("shrinkMatrix") }
)

setGeneric(name = "sensitivityAnalysis",
           def  = function(lp, ...) { standardGeneric("sensitivityAnalysis") }
)

setGeneric(name = "setColsNames",
           def  = function(lp, j, names) { standardGeneric("setColsNames") }
)

setGeneric(name = "setObjDir",
           def  = function(lp, lpdir) { standardGeneric("setObjDir") }
)

setGeneric(name = "setRhsZero",
           def  = function(lp) { standardGeneric("setRhsZero") }
)

setGeneric(name = "setRowsNames",
           def  = function(lp, i, names) { standardGeneric("setRowsNames") }
)

setGeneric(name = "setSolverParm",
           def  = function(lp, solverParm) { standardGeneric("setSolverParm") }
)

setGeneric(name = "singletonMetabolites",
           def  = function(object, ...) { standardGeneric("singletonMetabolites") }
)

setGeneric(name = "Snnz",
           def  = function(object) { standardGeneric("Snnz") }
)

setGeneric(name = "solveLp",
           def  = function(lp) { standardGeneric("solveLp") }
)

setGeneric(name = "solver",
           def  = function(object) { standardGeneric("solver") }
)
setGeneric(name = "solver<-",
           def  = function(object, value) { standardGeneric("solver<-") }
)

setGeneric(name = "status_code",
           def  = function(object) { standardGeneric("status_code") }
)
setGeneric(name = "status_code<-",
           def  = function(object, value) { standardGeneric("status_code<-") }
)

setGeneric(name = "status_meaning",
           def  = function(object) { standardGeneric("status_meaning") }
)
setGeneric(name = "status_meaning<-",
           def  = function(object, value) { standardGeneric("status_meaning<-") }
)

setGeneric(name = "status_num",
           def  = function(object) { standardGeneric("status_num") }
)
setGeneric(name = "status_num<-",
           def  = function(object, value) { standardGeneric("status_num<-") }
)

setGeneric(name = "subSys",
           def  = function(object) { standardGeneric("subSys") }
)
setGeneric(name = "subSys<-",
           def  = function(object, value) { standardGeneric("subSys<-") }
)

setGeneric(name = "uppbnd",
           def  = function(object) { standardGeneric("uppbnd") }
)
setGeneric(name = "uppbnd<-",
           def  = function(object, value) { standardGeneric("uppbnd<-") }
)

setGeneric(name = "uptake",
           def  = function(object) { standardGeneric("uptake") }
)
setGeneric(name = "uptake<-",
           def  = function(object, value) { standardGeneric("uptake<-") }
)

setGeneric(name = "uptReact",
           def  = function(object) { standardGeneric("uptReact") }
)

setGeneric(name = "uptMet",
           def  = function(object) { standardGeneric("uptMet") }
)

setGeneric(name = "verblevel",
           def  = function(object) { standardGeneric("verblevel") }
)
setGeneric(name = "verblevel<-",
           def  = function(object, value) { standardGeneric("verblevel<-") }
)

setGeneric(name = "writeProb",
           def  = function(lp, fname, ff = "lp", ...) { standardGeneric("writeProb") }
)


