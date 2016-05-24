#------------------------------------------------------------------------------#
#                             R interface to GLPK                              #
#------------------------------------------------------------------------------#

#  generics.R
#  R interface to GLPK.
#
#  Copyright (C) 2011-2014 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of glpkAPI.
#
#  GlpkAPI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  GlpkAPI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with glpkAPI  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                                   generics                                   #
#------------------------------------------------------------------------------#

setGeneric(name = "glpkPointer",
           def  = function(object) { standardGeneric("glpkPointer") }
)

setGeneric(name = "glpkPtrType",
           def  = function(object) { standardGeneric("glpkPtrType") }
)
setGeneric(name = "glpkPtrType<-",
           def  = function(object, value) { standardGeneric("glpkPtrType<-") }
)

setGeneric(name = "isNULLpointerGLPK",
           def  = function(object) { standardGeneric("isNULLpointerGLPK") }
)

setGeneric(name = "isGLPKpointer",
           def  = function(object) { standardGeneric("isGLPKpointer") }
)

setGeneric(name = "isTRWKSpointer",
           def  = function(object) { standardGeneric("isTRWKSpointer") }
)
