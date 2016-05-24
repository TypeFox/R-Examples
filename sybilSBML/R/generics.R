#------------------------------------------------------------------------------#
#                          Link to libSBML for sybil                           #
#------------------------------------------------------------------------------#

#  generics.R
#  Link to libSBML for sybil.
#
#  Copyright (C) 2010-2013 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybilSBML.
#
#  SybilSBML is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SybilSBML is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with SybilSBML.  If not, see <http://www.gnu.org/licenses/>.


#------------------------------------------------------------------------------#
#                                   generics                                   #
#------------------------------------------------------------------------------#

setGeneric(name = "sbmlPtrType",
           def  = function(object) { standardGeneric("sbmlPtrType") }
)

setGeneric(name = "sbmlPointer",
           def  = function(object) { standardGeneric("sbmlPointer") }
)

setGeneric(name = "sbmlDocKey",
           def  = function(object) { standardGeneric("sbmlDocKey") }
)

setGeneric(name = "sbmlFileName",
           def  = function(object) { standardGeneric("sbmlFileName") }
)

setGeneric(name = "isNULLpointerSBML",
           def  = function(object) { standardGeneric("isNULLpointerSBML") }
)

setGeneric(name = "isSBMLdocpointer",
           def  = function(object) { standardGeneric("isSBMLdocpointer") }
)

setGeneric(name = "isSBMLmodpointer",
           def  = function(object) { standardGeneric("isSBMLmodpointer") }
)

setGeneric(name = "sbmlInfos",
           def  = function(object) { standardGeneric("sbmlInfos") }
)

setGeneric(name = "sbmlWarnings",
           def  = function(object) { standardGeneric("sbmlWarnings") }
)

setGeneric(name = "sbmlErrors",
           def  = function(object) { standardGeneric("sbmlErrors") }
)

setGeneric(name = "sbmlFatals",
           def  = function(object) { standardGeneric("sbmlFatals") }
)

setGeneric(name = "getNumErrors",
           def  = function(object) { standardGeneric("getNumErrors") }
)

setGeneric(name = "printSlot",
           def  = function(object, ws) { standardGeneric("printSlot") }
)

