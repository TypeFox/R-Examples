#  optsol_robAnaClass.R
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


# optsol_robAnaClass


#------------------------------------------------------------------------------#
#                  definition of the class optsol_robAna                       #
#------------------------------------------------------------------------------#

setClass("optsol_robAna",
         representation(
              ctrlr   = "reactId",   # id of the control reaction,
              ctrlfl  = "numeric"     # fixed flux values of control reaction
         ),
         contains = "optsol_optimizeProb"
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# ctrlr
setMethod("ctrlr", signature(object = "optsol_robAna"),
          function(object) {
              return(object@ctrlr)
          }
)

setReplaceMethod("ctrlr", signature = (object = "optsol_robAna"),
                 function(object, value) {
                     object@ctrlr <- value
                     return(object)
                 }
)


# ctrlfl
setMethod("ctrlfl", signature(object = "optsol_robAna"),
          function(object) {
              return(object@ctrlfl)
          }
)

setReplaceMethod("ctrlfl", signature = (object = "optsol_robAna"),
                 function(object, value) {
                     object@ctrlfl <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

setMethod("plot", signature(x = "optsol_robAna", y = "missing"),
          function(x, y,
                   xlab = paste("Control Flux:", react_id(ctrlr(x))),
                   ylab = paste("Objective Function:", obj_func(x)),
                   type = "b",
                   pch = 20,
                   fillColorBg = "grey",
                   fillBg = TRUE,
                   absCtrl = TRUE,
                   ...) {

              if (isTRUE(absCtrl)) {
                  cr <- abs(x@ctrlfl)
              }
              else {
                  cr <- x@ctrlfl
              }

              plot(cr, x@lp_obj, type = "n", xlab = xlab, ylab = ylab)
              if (isTRUE(fillBg)) {
                  polygon(c(cr[1], cr, cr[length(cr)]),
                          c(min(x@lp_obj), x@lp_obj, min(x@lp_obj)),
                          col = fillColorBg, border = NA)
              }
              points(cr, x@lp_obj, type = type, pch = pch, ...)
          }
)

