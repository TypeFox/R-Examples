# Use multilinear interpolation on D-dimensional grid
# very inefficient, since grid is re-constructed for every interpolation request
# grid dimension (D) is inferred from length(location)
interpolateWithinGrid <-
function(	grid_string, 	# pre-computed string representing rectangular D-dimensional grid. Of format <list of x1 values>,<list of x2 values>,..,<list of xD values>,<list of z values in row-major format>   where each list is space-delimited and must be strictly increasing
			location, 		# new location at which to interpolate (vector of size D)
			default_z){		# value to use if requested location is outside of grid domain. If NULL, flat extrapolation is done (i.e. closest grid value is returned)
	# parse grid
	D 			= length(location);
	parts 		= strsplit(grid_string,split="[[:space:]]*,{1}[[:space:]]*",fixed=FALSE)[[1]];
	grid_points	= lapply(parts[1:D], FUN=function(p){ as.numeric(strsplit(p,split="[[:space:]]+",fixed=FALSE)[[1]]) } );
	z_values 	= as.numeric(strsplit(parts[D+1],split="[[:space:]]+",fixed=FALSE)[[1]]);
	N			= sapply(grid_points, length);
	if(length(z_values)!=prod(N)){
		stop(sprintf("ERROR parsing grid string: Grid size = %d = %s, but found %d z-values\n",prod(N),paste(as.character(N),collapse=' x '), length(z_values)));
	}
	EPSILON		= 1e-5; # relative grid tolerance, mainly to circumvent rounding errors
	
	if(any(N==0)){ return(default_z); } # no grid available. Returned default_z (even if NULL)
	location_du = location + abs(location)*EPSILON;
	location_dd = location - abs(location)*EPSILON;
	enclosing_box = sapply(c(1:D), FUN = function(d){ 
											if(location[d]<grid_points[[d]][1]){ return(0); }
											else if(location[d]>grid_points[[d]][N[d]]){ return(N[d]+1); }
											else{
												if(N[d]==1){
													return(0);
												}else{
													for(i in 2:(N[d])){ 
														if(grid_points[[d]][i-1]<=location_du[d] && grid_points[[d]][i]>=location_dd[d]){ return(i-1); } 
													}
													return(0); 
												}
											} 
										});
	if(!is.null(default_z)){
		if(any(enclosing_box<=0 | enclosing_box>N)) return(default_z);
	}
	box_volume_fractions = sapply(c(1:D), FUN = function(d){ 
													if(enclosing_box[d]<=0 | enclosing_box[d]>N[d]){ return(0); }
													else{ return((location[d]-grid_points[[d]][enclosing_box[d]])/(grid_points[[d]][enclosing_box[d]+1] - grid_points[[d]][enclosing_box[d]])); }
												});
	
	Nvertices = 2^D;
	z = 0;
	for(v in 0:(Nvertices-1)){
		vertex_index = 1; vertex_index_factor = 1;
		vertex_binary = as.numeric(intToBits(v)); vertex_binary = vertex_binary[1:D];
		vertex_weight = prod(sapply(1:D, FUN = function(d){ if(vertex_binary[d]==0){ return(1-box_volume_fractions[d]); }else{ return(box_volume_fractions[d]); } } ));
		for(d in D:1){
			if(enclosing_box[d]<=0){ vertex_index = vertex_index + 0; }
			else if(enclosing_box[d]>N[d]){ vertex_index = vertex_index + (N[d]-1)*vertex_index_factor; }
			else if(vertex_binary[d]==0){ vertex_index = vertex_index + (enclosing_box[d]-1)*vertex_index_factor; }
			else{ vertex_index = vertex_index + enclosing_box[d]*vertex_index_factor; }
			vertex_index_factor = vertex_index_factor * N[d];
		}
		z = z + z_values[vertex_index] * vertex_weight;
	}
	return(z);
}
