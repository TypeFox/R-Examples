inla.models=function(){
  if(requireNamespace("INLA", quietly=TRUE)){
    return(INLA::inla.models())
  } else {
    return(NULL)
  }
}