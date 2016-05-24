## this script defines a S4 class for Watershed objects
## writen by J.A. Torres-Matallana
## Institute for Geoinformatics, ifgi
## University of Muenster, Germany
## date: 9.08.2013

##  S4 class Watershed

setClass("Watershed", representation(station="SpatialPoints", subbasin="SpatialPolygonsDataFrame", zhyd="SpatialPolygonsDataFrame", river="SpatialLinesDataFrame", 	c1="SpatialPolygonsDataFrame", node="SpatialPointsDataFrame"))
		
