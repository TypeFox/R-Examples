f.get.marker.names <- function(data, n.vars){
##
## EXTRACT MARKER NAMES
##
#
if(class(data) == "gwaa.data"){
	## GWAA OBJECT
	.marker.names <- snpnames(data)
}else{
	## HAPLIN DATA MATRIX
	.marker.names <- colnames(data)
	if(n.vars > 0) .marker.names <- .marker.names[-(1:n.vars)]
	## REMOVE "l_" AT START AND "_m2" (AND SUCH) AT END
	.marker.names <- substring(.marker.names, first = 3, last = nchar(.marker.names) - 3)
	.marker.names <- unique(.marker.names)
}
#
return(.marker.names)
}
