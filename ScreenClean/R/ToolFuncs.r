




#NormalizeXY<-function(x,y) {
  ##################################################################
  #
  # This function normalized the data,
  # and calculate the inputs of the graphlet screening functions.
  #
  # Args:
  #   x: the predictor matrix, an n by p matrix
  #   y: the response, a length n vector
  #
  # Returns:
  #   x.mean: the means of each predictor
  #   y.mean: the mean of the response
  #   X: the normalized predictor matrix
  #   y.tilde: t(X)%*%y
  #   normalizer: the inverse of sqrt(p-1) times the stanard deviation of the predictors.
  #   gram: the normalized gram matrix.
  #
  ####################################################################
#  n<-length(y)
#  x.mean<-colMeans(x)
#  y.mean<-mean(y)
#  x<- x-matrix(1,n,1)%*%x.mean
#  y<-y-y.mean
#  gram.gram<-t(x) %*% x
#  normalizer.gram<-1/sqrt(diag(gram.gram))
#  gram<-diag(normalizer.gram)%*% gram.gram %*% diag(normalizer.gram)
#  x.n<-x %*% diag(normalizer.gram)
#  y.tilde<-t(x.n)%*%y
#  return(list(x.mean=x.mean,y.mean=y.mean,X=x.n,y.tilde=y.tilde,normalizer=normalizer.gram, gram=gram))
#}


ThresholdGram<- function(gram.full,delta=1/log(dim(gram.full)[1])) {
  #####################################################################
  #
  # This is a simple function that threshold the gram matrix
  # 
  # Args:
  #   gram.full: the original gram matrix
  # 
  # Returns:
  #   gram: the gram matrix after thresholding
  #   gram.bias: the bias that the thresholding makes
  #
  #######################################################################
 
  gram.bias<-gram.full*(abs(gram.full)<delta)
  gram.sd<-gram.full-gram.bias
  return(list(gram=gram.sd,gram.bias=gram.bias))
}





