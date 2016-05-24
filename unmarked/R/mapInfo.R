

setClass("mapInfo",
		representation(coordinates = "matrix", 
				projection = "optionalCharacter", 
				parameters = "optionalNumeric", 
				orientation = "optionalNumeric"))

# note order: long, lat

mapInfo <- function(coordinates, projection = NULL, parameters = NULL, orientation = NULL) {
	new("mapInfo", coordinates = coordinates, projection = projection, parameters = parameters, orientation = orientation)
}


setClassUnion("optionalMapInfo", c("mapInfo","NULL"))