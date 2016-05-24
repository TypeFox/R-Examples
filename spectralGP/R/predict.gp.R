"predict.gp" <-
function(object,newdata=NULL,mapping=NULL,...){
  # return the process values on the grid or at the gridpoints nearest to a set of supplied locations or a supplied mapping from locations to gridpoints 
  if(is.null(newdata) & is.null(mapping)){
    # return process values on grid
    m1=object$gridsize[1]
    m2=object$gridsize[2]
    if(object$d==2){
      return(matrix(object$process,m1,m2)[1:(m1/2+1),1:(m2/2+1)])
    } else{
      return(matrix(object$process,m1,m2)[1:(m1/2+1),])
    }
  }
  if(!is.null(mapping) & !is.null(newdata)){
    warning("Both newdata and mapping are supplied; newdata are being ignored")
  }
  if(is.null(mapping)){
    # need to calculate mapping before doing prediction
    mapping=new.mapping(object,newdata)
  }
  if(min(mapping<1) || max(mapping)>prod(object$gridsize)){
    stop(" 'mapping' must be a vector that legitimately maps process locations (on a grid) to a set of locations")
  }
  return(object$process[mapping])
}
