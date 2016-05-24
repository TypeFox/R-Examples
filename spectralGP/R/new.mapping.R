"new.mapping" <-
function(object,locations){
  # finds the nearest grid point to each location in locations and puts that in a new component, 'mapping'
  # this assumes locations are in (0,1)^d
  # returns vector that maps locations to grid locations
  if(object$d==1){
    if(!
       (is.vector(locations) || (is.array(locations) && length(dim(locations))==1) || (is.matrix(locations) && min(dim(locations))==1)
        )){
      stop("Locations must be one-dimensional")
    }
  } else{
    if(length(dim(locations))!=2 || nrow(locations)<1 || ncol(locations)!=2){
      stop("locations must be a two-column matrix-like object, where the first column indicates the x-coordinate and the second the y-coordinate")
    }
  }
  if(min(locations)<0 || max(locations)>1){
    warning("locations should lie in or (for test locations) near (0,1)^2.\n")
  }
  if(object$d==1){
    mapping=round((object$gridsize[1]/2)*locations)+1
  } else{
    xdim=round((object$gridsize[1]/2)*locations[,1])+1
    ydim=round((object$gridsize[2]/2)*locations[,2])+1
    mapping=(ydim-1)*object$gridsize[2]+xdim  
  }
  return(mapping)
}
