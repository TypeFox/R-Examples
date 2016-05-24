#  optObj_pointer.R
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
#             class definitions for pointers to problem objects                #
#------------------------------------------------------------------------------#

setClass(Class = "lpExtPtr", contains = "VIRTUAL")     # lpSolveAPI
setClass(Class = "glpkPtr",  contains = "VIRTUAL")     # glpkAPI
setClass(Class = "clpPtr",   contains = "VIRTUAL")     # clpAPI
setClass(Class = "cplexPtr", contains = "VIRTUAL")     # cplexAPI


#------------------------------------------------------------------------------#
#                     definition of the class cplexPointer                     #
#------------------------------------------------------------------------------#

setClass(Class = "cplexPointer",
         representation(
             env = "cplexPtr",
             lp  = "cplexPtr"
         ),
)


#------------------------------------------------------------------------------#
#                            default constructor                               #
#------------------------------------------------------------------------------#

# contructor for class cplexPointer
setMethod(f = "initialize",
          signature = "cplexPointer",
          definition = function(.Object, en, pr) {

              if ( (!missing(en)) || (!missing(pr)) ) {
                  if ( (cplexAPI::isCPLEXenvPointer(en)) &&
                       (cplexAPI::isCPLEXprobPointer(pr)) ) {
                  
                      .Object@env <- en
                      .Object@lp  <- pr

                  }
              }
              return(.Object)
          }
)


#------------------------------------------------------------------------------#
#                    definition of the class pointerToProb                     #
#------------------------------------------------------------------------------#

# pointer representation in class optObj
setClassUnion(name    = "pointerToProb",
              members = c("externalptr",
                          "lpExtPtr",
                          "glpkPtr",
                          "clpPtr",
                          "cplexPointer"
             )
)


