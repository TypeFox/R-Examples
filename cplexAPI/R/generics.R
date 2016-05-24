#------------------------------------------------------------------------------#
#                     R Interface to C API of IBM ILOG CPLEX                   #
#------------------------------------------------------------------------------#

#  generics.R
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
#                                   generics                                   #
#------------------------------------------------------------------------------#

setGeneric("summary")


setGeneric(name = "err",
           def  = function(object) { standardGeneric("err") }
)

setGeneric(name = "errmsg",
           def  = function(object) { standardGeneric("errmsg") }
)

setGeneric(name = "errnum",
           def  = function(object) { standardGeneric("errnum") }
)
setGeneric(name = "errnum<-",
           def  = function(object, value) { standardGeneric("errnum<-") }
)

setGeneric(name = "cplexPointer",
           def  = function(object) { standardGeneric("cplexPointer") }
)

setGeneric(name = "cplexPtrType",
           def  = function(object) { standardGeneric("cplexPtrType") }
)
setGeneric(name = "cplexPtrType<-",
           def  = function(object, value) { standardGeneric("cplexPtrType<-") }
)

setGeneric(name = "isCPLEXprobPointer",
           def  = function(object) { standardGeneric("isCPLEXprobPointer") }
)

setGeneric(name = "isCPLEXenvPointer",
           def  = function(object) { standardGeneric("isCPLEXenvPointer") }
)

setGeneric(name = "isCPLEXfilePointer",
           def  = function(object) { standardGeneric("isCPLEXfilePointer") }
)

setGeneric(name = "isCPLEXchanPointer",
           def  = function(object) { standardGeneric("isCPLEXchanPointer") }
)

setGeneric(name = "isCPLEXtermPointer",
           def  = function(object) { standardGeneric("isCPLEXtermPointer") }
)

setGeneric(name = "isNULLpointerCPLEX",
           def  = function(object) { standardGeneric("isNULLpointerCPLEX") }
)
