simprepSpatial <-
function(regions, coord, radius=NULL, form=c("cube", "sphere","manual"), fading=0){

	
	if(missing(form)){
		form <- "cube"
	}
	if(!is.list(coord)&&(regions>1)){
		stop("Coord should be a list. See ?simprepSpatial for further details.")
	}else{
	  if(is.numeric(coord)){
		coord <- list(coord)
	  }
	}
	if(length(coord)!=regions){
		stop("Mismatch between coord and regions.  See ?simprepSpatial for further details.")
	}
	if(is.null(radius)&&(form!="manual")){
		stop("Argument radius is not specified. See ?simprepSpatial for further details.")
	}
	if(length(radius)==1){
		radius <- rep(radius, regions)
	} else if(length(radius)!=regions){
		stop("Mismatch between radius and regions. See ?simprepSpatial for further details.")
	}
	if(length(form)==1){
		form <- rep(form, regions)
	}
	if(length(fading)==1){
		fading <- rep(fading, regions)
	}
	
	spat <- list()
	for(i in 1:regions){
		spat[[i]] <- list()
		name <- paste("region",i, sep="")
		spat[[i]]$name <- name
		spat[[i]]$coord <- coord[[i]]
		spat[[i]]$radius <- radius[i]
		spat[[i]]$form <- form[i]
		spat[[i]]$fading <- fading[i]
	}
	return(spat)
}

