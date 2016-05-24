# An equivalent to the Matlab function of the same name, for internal use
# @param nrows Number of rows (default: 1)
# @param ncols Number of columns (default: nrows)
# @return A matrix of zeros with nrows rows and ncols columns
zeros = function(nrows=1,ncols=nrows) {
  return(matrix(0,nrows,ncols));
}

# An equivalent to the Matlab function of the same name, for internal use
# @param nrows Number of rows (default: 1)
# @param ncols Number of columns (default: nrows)
# @return A matrix of ones with nrows rows and ncols columns
ones = function(nrows=1,ncols=nrows) {
  return(matrix(1,nrows,ncols));
}

# Sigmoid function, for internal use
# @param x
# @return sigmoid(x)
sigmoid = function(x) {
  return(1/(1 + exp(-x)));
}

# Gradient sigmoid function, for internal use
# @param x
# @return Gradient for gradient descent on sigmoid(x)
sigmoidGradient = function(z) {
  return((1-1/(1+exp(-z)))/(1+exp(-z)));
}

# Random initilization of weights, for internal use
# @param L_in Number of input variables
# @param L_out Number of output variables
# @return Gradient for gradient descent on sigmoid(x)
randInitializeWeights = function(L_in, L_out) {
  epsilon_init = 0.12;
  matrix(runif(L_out*(L_in+1)),L_out,L_in+1) * 2 * epsilon_init - epsilon_init;
}


# normalize a vector to a [0;1] range
normalize=function(obj) {
  denom=max(obj,na.rm=T)-min(obj,na.rm=T);
  if(denom==0) return(rep(0,length(obj)));
  return((obj-min(obj,na.rm=T))/denom);
}
#' Normalize data
#' @description Normalize all columns of a dataframe so that all values are in [0;1] and for each column the maximum
#' value is 1 and the minimum 0.
#' 
#' newx=(x-min(X))/(max(X)-min(X))
#' @param dframe The dataframe to be normalized
#' @return The normalized dataframe
#' @export
sN.normalizeDF=function(dframe) {
  nVar=ncol(dframe);
  out=dframe;
  for(i in 1:nVar){
    out[[i]]=normalize(out[[i]]);
  }
  return(out);
}