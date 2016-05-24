# TODO: Add comment
# 
# Author: cws
###############################################################################

attribute_container_as_list <- function(attributeContainer){
	.Call( "attribute_container_as_list", attributeContainer, PACKAGE = "NetSim" )
}

as.numeric.NetSimAttributeContainer <- function(x, ...){
	return(attribute_container_as_list(x))
}

as.double.NetSimAttributeContainer <- function(x, ...){
	return(attribute_container_as_list(x))
}


create_attribute_container <- function(numericVector){
	.Call( "create_attribute_container", numericVector, PACKAGE = "NetSim" )
}

create_scale_attribute_container <- function(numericVector, min = 0, max = 1, by = 1){
	.Call("create_scale_attribute_container", 
			numericVector, min, max, by, PACKAGE = "NetSim")
}

set_value <- function(attributeContainer, i, value){
	.Call("set_value", attributeContainer, i, value, PACKAGE = "NetSim")
}

print.NetSimAttributeContainer <- function(x, ...){
	print("Pointer to NetSim attribute container object:", quote = FALSE)
	print(attribute_container_as_list(x))
}

