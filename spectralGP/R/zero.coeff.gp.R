"zero.coeff.gp" <-
function(object,...){
  # zeroes out the coefficients
  object$coeff=matrix(0,nrow=object$gridsize[1],ncol=object$gridsize[2])
  updateprocess(object)
  return(NULL)
}
