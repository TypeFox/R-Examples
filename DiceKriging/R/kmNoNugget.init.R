`kmNoNugget.init` <-
function(model, fn, fnscale) {

	parinit <- model@parinit
  ninit <- model@control$pop.size
	param.n <- model@covariance@param.n      
	
	if (length(parinit) > 0) {
	  matrixinit <- matrix(parinit, nrow = param.n, ncol = 1)    
	} else {
	  lower <- model@lower
	  upper <- model@upper
	  if (existsMethod("paramSample", signature = class(model@covariance))) {
	    matrixinit <- paramSample(model@covariance, n=ninit, lower=lower, upper=upper, y=model@y)
	  } else {
	    # sample ninit design points, generated from uniform [lower, upper]
	    matrixinit <- matrix(runif(ninit*param.n), nrow = param.n, ncol = ninit)
	    matrixinit <- lower + matrixinit*(upper - lower)
	  }
	}
  
  # take the best point(s)
	fninit <- apply(matrixinit, 2, fn, model)
  selection <- sort(fninit, decreasing = (fnscale < 0), index.return = TRUE)$ix
  selection <- selection[1:model@control$multistart]
	parinit <- matrixinit[, selection, drop = FALSE]  
  # for one point : parinit <- matrixinit[, which.max(fninit), drop = FALSE]  		
	
	covinit <- list()
  for (i in 1:model@control$multistart){
	  covinit[[i]] <- vect2covparam(model@covariance, parinit[,i])
	}
  
  return(list(par = parinit, 
              value = fninit[selection], 
              cov = covinit,
              lower = model@lower,
              upper = model@upper))	       
	
}	
