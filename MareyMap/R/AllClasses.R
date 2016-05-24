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


#--- Interpolation class -------------------------------------------------------

setClass("Interpolation",
#  contains = "VIRTUAL",
	representation(
    name = "character",
		color = "character",
		physicalPositions = "vector",
		rates = "vector",
		visible = "logical",
		persistent = "logical"
	),
	prototype(
    name = "name",
		color = "#969696",
		visible = T,
		rates = vector(),
		physicalPositions = vector(),
		persistent = T
	)
)


#--- InterpolationParamClass ---------------------------------------------------

setClass("InterpolationParam",
	representation(
		paramName = "character",
		paramType = "character",
		paramDesc = "character",
		paramDefault = "ANY",
		
		paramValues = "vector",
		paramMin = "ANY",
		paramMax = "ANY",
		paramFun = "character"
		),
	prototype(
		paramName = "parameter's name",
		paramType = "parameter's type",
		paramDesc = "a short description of the parameter",
		paramDefault = "ANY",
		paramValues = NULL,
		paramMin = NULL,
		paramMax = NULL,
		paramFun = NULL)
)


#--- MareyMap ------------------------------------------------------------------

setClass("MareyMap", 
	representation(
		setName = "character", 
		mapName = "character",  
		markerNames = "vector",
		physicalPositions = "vector", 
		geneticDistances = "vector",
		markerValidity = "vector",
		interpolations = "list"
	)
)


#--- MapSet class ---------------------------------------------------------

setClass("MapSet",
	representation(
		maps = "list",
		setName = "character"
	)
)


#--- MapCollection -------------------------------------------------------------

setClass("MapCollection",
	representation(
		sets = "list"
	)
)
