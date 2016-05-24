"change.param.gp" <-
function(object,new.specdens.param=NULL,new.variance.param=NULL,...){
  if(!is.null(new.specdens.param) & length(new.specdens.param)!=length(object$specdens.param)){
    stop('New spectral density parameter vector is not of the proper length: ',length(object$specdens.param),'\n')
  }
  if(!is.null(new.variance.param) & length(new.variance.param)!=1){
    stop('New variance parameter should be a scalar.\n')
  } 
  if(!is.null(new.specdens.param)){
    object$specdens.param=new.specdens.param
  }
  if(!is.null(new.variance.param)){
    object$variance.param=new.variance.param
  }
  calc.variances(object)
  return(NULL)
}
