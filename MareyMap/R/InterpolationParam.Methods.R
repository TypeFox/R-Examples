# Copyright 2006 Laboratoire de Biologie et de Biometrie Appliquée 
# (UMR 5558);CNRS; Univ. Lyon 1, 43 bd 11 nov, 69622, 
# Villeurbanne Cedex, France.
#
# This file is part of MareyMap.
#
# MareyMap is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# MareyMap is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MareyMap; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


#--- constructor
setMethod("InterpolationParam", "missing", function(dummy) {
	new("InterpolationParam")
})


#--- accessors
setMethod("paramName", "InterpolationParam", function(object) object@paramName)
setMethod("paramType", "InterpolationParam", function(object) object@paramType)
setMethod("paramDesc", "InterpolationParam", function(object) object@paramDesc)	
setMethod("paramValues", "InterpolationParam", function(object) object@paramValues)
setMethod("paramMin", "InterpolationParam", function(object) object@paramMin)
setMethod("paramMax", "InterpolationParam", function(object) object@paramMax)
setMethod("paramDefault", "InterpolationParam", function(object) object@paramDefault)
setMethod("paramFun", "InterpolationParam", function(object) object@paramFun)	


#--- replace methods
setReplaceMethod("paramName",	c("InterpolationParam", "character"), function(object, value) {
	object@paramName <- value
	object
})

setReplaceMethod("paramType",	c("InterpolationParam", "character"), function(object, value) {
	object@paramType <- value
	object
})

setReplaceMethod("paramDesc",	c("InterpolationParam", "character"),	function(object, value) {
	object@paramDesc <- value
	object
})	

setReplaceMethod("paramValues", c("InterpolationParam", "vector"), function(object, value) {
	object@paramValues <- value
	object
})

setReplaceMethod("paramMin", c("InterpolationParam", "ANY"), function(object, value) {
	object@paramMin <- value
	object
})

setReplaceMethod("paramMax", c("InterpolationParam", "ANY"), function(object, value) {
	object@paramMax <- value
	object
})

setReplaceMethod("paramDefault", c("InterpolationParam", "ANY"), function(object, value) {
	object@paramDefault <- value
	object
})

setReplaceMethod("paramFun", c("InterpolationParam", "character"), function(object, value) {
	object@paramFun <- value
	object
})

