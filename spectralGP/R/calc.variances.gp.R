"calc.variances.gp" <-
function(object,...){
  # calculate prior variances based on specdens.parameters; note that if only variance.param changes, this unnecessarily recalculates spectral density
  variances=object$variance.param*matrix(object$specdens(object$omega,object$specdens.param,d=object$d),nrow=object$gridsize[1],ncol=object$gridsize[2],byrow=FALSE)
  variances=0.5*variances*prod(object$gridsize) # ensures that variances are on the correct scale for doing Gibbs sampling via the Wikle (2002) algorithm
  # next lines set variances for real coefficients back to value from spectral density
  variances[1,1]=2*variances[1,1]
  variances[(object$gridsize[1]/2+1),1]=2*variances[(object$gridsize[1]/2+1),1]
  if(object$d==2){
    variances[1,(object$gridsize[2]/2+1)]=2*variances[1,(object$gridsize[2]/2+1)]
    variances[(object$gridsize[1]/2+1),(object$gridsize[2]/2+1)]=2*variances[(object$gridsize[1]/2+1),(object$gridsize[2]/2+1)]
  }
  object$variances=variances
  return(NULL)
}
