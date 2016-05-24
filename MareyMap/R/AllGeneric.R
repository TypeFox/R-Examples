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


setGeneric("argList", function(object) standardGeneric("argList"))
setGeneric("color", function(object) standardGeneric("color"))
setGeneric("color<-", function(object, value) standardGeneric("color<-"))
setGeneric("createOrder", function(object) standardGeneric("createOrder"))
setGeneric("geneticDistances", function(object) standardGeneric("geneticDistances"))
setGeneric("geneticDistances<-", function(object, value) standardGeneric("geneticDistances<-"))
setGeneric("getInterList", function(dummy) standardGeneric("getInterList"))
setGeneric("interpolate", function(object, map) standardGeneric("interpolate"))
setGeneric("interpolation",	function(object, inter_name) standardGeneric("interpolation"))
setGeneric("interpolation<-", function(object, inter_name, value) standardGeneric("interpolation<-"))
setGeneric("InterpolationParam", function(dummy) standardGeneric("InterpolationParam"))
setGeneric("interpolations", function(object) standardGeneric("interpolations"))
setGeneric("interpolations<-", function(object, value) standardGeneric("interpolations<-"))
setGeneric("MapCollection",function(x) standardGeneric("MapCollection"))
setGeneric("mapName", function(object) standardGeneric("mapName"))
setGeneric("mapName<-", function(object, value)	standardGeneric("mapName<-"))
setGeneric("mapNames", function(object) standardGeneric("mapNames"))
setGeneric("MareyMap", function(data_table, column_names = list(name = "name", phys = "phys", gen = "gen"),	set_name = "", map_name = "")	standardGeneric("MareyMap"))
setGeneric("markerNames",	function(object) standardGeneric("markerNames"))
setGeneric("markerNames<-",	function(object, value)	standardGeneric("markerNames<-"))
setGeneric("markerValidity",function(object) standardGeneric("markerValidity"))
setGeneric("markerValidity<-", function(object, value) standardGeneric("markerValidity<-"))
setGeneric("name", function(object) standardGeneric("name"))
setGeneric("name<-", function(object,value) standardGeneric("name<-"))
setGeneric("paramDefault",function(object) standardGeneric("paramDefault"))
setGeneric("paramDefault<-", function(object, value) standardGeneric("paramDefault<-"))
setGeneric("paramDesc", function(object) standardGeneric("paramDesc"))	
setGeneric("paramDesc<-",	function(object, value) standardGeneric("paramDesc<-"))
setGeneric("paramFun", function(object) standardGeneric("paramFun"))	
setGeneric("paramFun<-", function(object, value) standardGeneric("paramFun<-"))
setGeneric("paramMax", function(object) standardGeneric("paramMax"))
setGeneric("paramMax<-", function(object, value) standardGeneric("paramMax<-"))
setGeneric("paramMin", function(object) standardGeneric("paramMin"))
setGeneric("paramMin<-", function(object, value) standardGeneric("paramMin<-"))
setGeneric("paramName", function(object) standardGeneric("paramName"))
setGeneric("paramName<-", function(object, value) standardGeneric("paramName<-"))
setGeneric("paramType", function(object) standardGeneric("paramType"))
setGeneric("paramType<-", function(object, value) standardGeneric("paramType<-"))
setGeneric("paramValues", function(object) standardGeneric("paramValues"))
setGeneric("paramValues<-", function(object, value) standardGeneric("paramValues<-"))
setGeneric("persistent", function(object) standardGeneric("persistent"))
setGeneric("persistent<-", function(object, value) standardGeneric("persistent<-"))
setGeneric("physicalPositions", function(object) standardGeneric("physicalPositions"))
setGeneric("physicalPositions<-", function(object, value)	standardGeneric("physicalPositions<-"))
setGeneric("plotMarkers", function(object, ...) standardGeneric("plotMarkers"))
setGeneric("plotModel", function(object, ...) standardGeneric("plotModel"))
setGeneric("plotModels", function(object, ...) standardGeneric("plotModels"))
setGeneric("plotRate", function(object, ...) standardGeneric("plotRate"))
setGeneric("plotRates", function(object, ...) standardGeneric("plotRates"))
setGeneric("query", function(object, pos) standardGeneric("query"))
setGeneric("rates", function(object) standardGeneric("rates"))
setGeneric("rates<-", function(object, value) standardGeneric("rates<-"))
setGeneric("registerInterpolationMethod", function(name, classname) standardGeneric("registerInterpolationMethod"))
setGeneric("removeMarker", function(object, value) standardGeneric("removeMarker"))
setGeneric("MapSet", function(set_name) standardGeneric("MapSet"))
setGeneric("setName", function(object) standardGeneric("setName"))
setGeneric("setName<-", function(object, value)	standardGeneric("setName<-"))
setGeneric("setNames", function(object) standardGeneric("setNames"))
setGeneric("textFile", function(object,file) standardGeneric("textFile"))
setGeneric("userParam", function(object) standardGeneric("userParam"))
setGeneric("validPositions", function(object) standardGeneric("validPositions"))
setGeneric("visible", function(object) standardGeneric("visible"))
setGeneric("visible<-", function(object,value) standardGeneric("visible<-"))
