#  optsol_genedelClass.R
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


# optsol_genedelClass


#------------------------------------------------------------------------------#
#                  definition of the class optsol_genedel                      #
#------------------------------------------------------------------------------#

setClass("optsol_genedel",
           representation(
                fluxdels  = "list",
                hasEffect = "logical"
           ),
           contains = "optsol_fluxdel"
)


#------------------------------------------------------------------------------#
#                            setters and getters                               #
#------------------------------------------------------------------------------#

# reactions
setMethod("fluxdels", signature(object = "optsol_genedel"),
          function(object) {
              return(object@fluxdels)
          }
)

setReplaceMethod("fluxdels", signature = (object = "optsol_genedel"),
                 function(object, value) {
                     object@fluxdels <- value
                     return(object)
                 }
)


# hasEffect
setMethod("hasEffect", signature(object = "optsol_genedel"),
          function(object) {
              return(object@hasEffect)
          }
)

setReplaceMethod("hasEffect", signature = (object = "optsol_genedel"),
                 function(object, value) {
                     object@hasEffect <- value
                     return(object)
                 }
)


#------------------------------------------------------------------------------#
#                               other methods                                  #
#------------------------------------------------------------------------------#

#setMethod("ind2id", signature = (object = "optsol_genedel"),
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
#                         "fluxdels" = {
#                             out <- lapply(fluxdels(object),
#                                           function(x) react_id(object)[x]
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


setMethod("deleted", signature = (object = "optsol_genedel"),
                 function(object, i) {
                     value <- fluxdels(object)[[i]]
                     return(value)
                 }
)


