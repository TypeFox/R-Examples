ratioaaddon <-
function(params, x, batch) {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.")

  if(class(params) != "ratioa")
     stop("Input parameter 'params' has to be of class 'ratioa'.")
	 
  if(ncol(params$xadj) != ncol(x))
    stop("Number of variables in test data matrix different to that of training data matrix.")	 

  batches = levels(batch)
  nbatches = length(batches)

  means = as.list(rep(0,nbatches))
  xadj = x  
  for (i in 1:nbatches) {
    means[[i]] <- colMeans(x[batch==batches[i],])
	meanabovezero <- which(means[[i]]!=0)
    xadj[batch==batches[i],meanabovezero] = scale(x[batch==batches[i],meanabovezero],center=rep(0,ncol(x[,meanabovezero])),scale=means[[i]][meanabovezero])
  }
  
  return(xadj)

}
