select.eigentraits <-
function(data.obj, traits.which = c(1,2)){
	
	ET <- data.obj$ET
	
	if(max(traits.which) > dim(ET)[2]){
		stop("Undefined ET have been specified.")
		}
	if(length(unique(traits.which)) < length(traits.which)){
		stop("There are redundant ET specified.")
		}
	selected.ET <- ET[,traits.which]
	data.obj$ET <- selected.ET
	return(data.obj)
	
}
