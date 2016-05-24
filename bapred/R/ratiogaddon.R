ratiogaddon <-
function(params, x, batch) {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!(is.factor(batch) & all(levels(batch)==(1:length(levels(batch))))))
    stop("'batch' has to be of class 'factor' with levels '1','2',...")  
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 

  if(class(params) != "ratiog")
    stop("Input parameter 'params' has to be of class 'ratiog'.")
	
  if(ncol(params$xadj) != ncol(x))
    stop("Number of variables in test data matrix different to that of training data matrix.")	 	

  batches = levels(batch)
  nbatches = length(batches)

  means = as.list(rep(0,nbatches))
  xadj = x  
  for (i in 1:nbatches) {
    xadj[batch==batches[i],] <- apply(x[batch==batches[i],], 2, function(x) {
	  if(all(x <= 0))
	    x <- rep(min(-(x[x < 0])), length(x))
	  else
    	    x[x<=0] <- min(x[x>0])
	  x
    })
    means[[i]] <- apply(xadj[batch==batches[i],], 2, function(x) exp(mean(log(x))))
    xadj[batch==batches[i],] = scale(xadj[batch==batches[i],],center=rep(0,ncol(xadj)),scale=means[[i]])
  }
  
  return(xadj)

}
