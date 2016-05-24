get.marker.idx <-
function(data.obj, markers){
	
	is.char <- as.logical(is.na(suppressWarnings(as.numeric(markers[1]))))
		
	if(is.char){
		marker.loc <- match(markers, data.obj$marker.names)
		return(marker.loc)
		}
	
	if(!is.char){
		marker.loc <- match(markers, data.obj$marker.num)
		return(marker.loc)
		}
	
	
}
