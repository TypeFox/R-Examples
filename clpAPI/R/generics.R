#------------------------------------------------------------------------------#
#                           R interface to COIN-OR Clp                         #
#------------------------------------------------------------------------------#

#  generics.R
#  R interface to COIN-OR Clp.
#
#  Copyright (C) 2011-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of clpAPI.
#
#  ClpAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ClpAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with clpAPI  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                                   generics                                   #
#------------------------------------------------------------------------------#

setGeneric(name = "clpPointer",
           def  = function(object) { standardGeneric("clpPointer") }
)

setGeneric(name = "clpPtrType",
           def  = function(object) { standardGeneric("clpPtrType") }
)
setGeneric(name = "clpPtrType<-",
           def  = function(object, value) { standardGeneric("clpPtrType<-") }
)

setGeneric(name = "isNULLpointerCLP",
           def  = function(object) { standardGeneric("isNULLpointerCLP") }
)

setGeneric(name = "isCLPpointer",
           def  = function(object) { standardGeneric("isCLPpointer") }
)

