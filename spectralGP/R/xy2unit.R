"xy2unit" <-
function(locations,locations.scale=NULL){
  # scales locations to (0,1)^d based on the range of values in the locations.scale matrix (or the locations matrix if no locations.scale matrix is supplied
  d=2
  if(is.vector(locations) || (is.array(locations) && length(dim(locations))==1) || (is.matrix(locations) && min(dim(locations))==1)){
    print("Scaling one-dimensional locations to (0,1)")
    d=1
  } else{
    if(length(dim(locations))!=2 || nrow(locations)<1 || ncol(locations)!=2){
      stop("Locations must be a vector or a matrix with two columns, where the first column indicates the x-coordinate and the second the y-coordinate")
    } else{
      print("Scaling two-dimensional locations to (0,1)^2")
    }
  }
  if(is.null(locations.scale)){
    locations.scale=locations
  }
  if(length(dim(locations.scale))!=length(dim(locations))){
    stop("Dimension of locations and locations.scale must be the same")
  }
  if(d==1){
    locations.scale=matrix(locations.scale,ncol=1)
    locations=matrix(locations,ncol=1)
  }
  maxs=apply(locations.scale,2,max)
  mins=apply(locations.scale,2,min)
  diffs=maxs-mins
  maxdiff=max(diffs)
  for(i in 1:d){
    locations[,i]=(locations[,i]-mins[i])/maxdiff
  }
  if(d==1){
    locations=c(locations)
  }
  return(locations)
}
