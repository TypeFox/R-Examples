svabaaddon <-
function(params, x) {

  if(any(is.na(x)))
	stop("Data contains missing values.") 
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 

  if(class(params) != "svatrain")
     stop("Input parameter 'params' has to be of class 'svatrain'.")
	 
  if(ncol(params$xadj) != ncol(x))
    stop("Number of variables in test data matrix different to that of training data matrix.")	 

  ##require("sva")

  if (!is.null(params$svobj)) {
    fsvaobj <- sva::fsva(dbdat = t(params$xtrain), mod = cbind(1, 
        as.numeric(params$ytrain) - 1), sv = params$svobj, newdat = t(x), 
        method = params$algorithm)
    xadj <- t(fsvaobj$new)
  }
  else
    xadj <- x

  return(xadj)

}
