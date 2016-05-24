################################################################################
# File: SpatialDesign-class
# Purpose: Create the definition for class SpatialDesign
# Programmer: Tom Kincaid
# Date: May 1, 2015
################################################################################

setClass(Class = "SpatialDesign",
	contains = "SpatialPointsDataFrame", 
	slots = c(design = "list"),
	validity = function(object) {
		if(!is.list(object@design))
			stop("The design argument to the SpatialDesign function must be a list.")
		if(length(names(object@design)) == 0)
			stop("The design argument to the SpatialDesign function must be named.")
		return(TRUE)
	}
)

SpatialDesign <- function(design, sp_obj) {
	new("SpatialDesign", design = design, data = sp_obj@data,
		coords.nrs = sp_obj@coords.nrs, coords = sp_obj@coords, bbox = sp_obj@bbox,
		proj4string = sp_obj@proj4string)
}
