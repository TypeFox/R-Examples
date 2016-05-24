nobaaddon <-
function(params, x, batch) {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 

  if(class(params) != "noba")
    stop("Input parameter 'params' has to be of class 'noba'.")
	
  if(ncol(params$xadj) != ncol(x))
    stop("Number of variables in test data matrix different to that of training data matrix.")	 	

  xadj <- x
	
  return(xadj)

}
