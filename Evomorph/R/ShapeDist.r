ShapeDist <- function(shapes, reference){
  if ((!(is.list(shapes)))&&(!(is.array(shapes))))
  {stop("Shapes must be list type or array type")}
  if (is.array(shapes))
  {
    shapesaux<-list();
    dimen<-dim(shapes)
    for (j in 1:(dimen[3])){
      shapesaux[[j]]=shapes[,,j];
    }
    shapes=shapesaux
  }
  if (length(shapes[[1]])!=length(reference))
  {stop("Shapes and reference must be the same length")}
  distances=array()
  for(i in 1:length(shapes)){distances[[i]]=sqrt(sum((shapes[[i]]-reference)^2))}
  return(distances)
}