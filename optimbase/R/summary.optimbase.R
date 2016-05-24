# Copyright (C) 2010-2014 - Sebastien Bihorel
#
# This file must be used under the terms of the CeCILL.
# This source file is licensed as described in the file COPYING, which
# you should have received as part of this distribution. The terms
# are also available at
# http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
#

summary.optimbase <- function(object,showhistory,...){

  if(missing('showhistory'))
    showhistory <- FALSE
  
  # 
  # Parameters estimates and boundaries
  #
  
  cat(sprintf('\nNumber of Estimated Variable(s): %d\n', object$numberofvariables))
  
  if (object$numberofvariables>0){
    estimates <- data.frame(object$x0,object$xopt)
    names(estimates) <- c('Initial','Final')
    if (length(object$boundsmin)>0 && !is.null(object$boundsmin)){
      tmp <- names(estimates)
      estimates <- cbind(estimates,object$boundsmin)
      names(estimates) <- c(tmp, "Lower bound")
    }
    if (length(object$boundsmax)>0 && !is.null(object$boundsmax)){
      tmp <- names(estimates)
      estimates <- cbind(estimates,object$boundsmax)
      names(estimates) <- c(tmp, "Upper bound") 
    }
    cat("\nEstimated Variable(s):\n")
    print(estimates)
  }
  
  if (object$nbineqconst>=1){
    cat(sprintf('\nNumber of Inequality Contraints: %s\n', object$nbineqconst))
  }
  
  #
  #Optimization info
  #
  cat('\nCost Function:\n')
  print(object$fun)
  cat('\nCost Function Argument(s):\n')
  print(object$costfargument)
  
  cat('\nOptimization:\n') 
  cat(sprintf('- Status: \"%s\"\n', object$status))
  cat(sprintf('- Initial Cost Function Value: %f\n', object$fx0))
  cat(sprintf('- Final Cost Function Value: %f\n', object$fopt))
  cat(sprintf('- Number of Iterations (max): %d (%d)\n', object$iterations,object$maxiter))
  cat(sprintf('- Number of Function Evaluations (max): %d (%d)\n\n', object$funevals,object$maxfunevals))
  
  if (object$storehistory && showhistory==TRUE){
    tmp <- as.data.frame(transpose(sapply(object$historyxopt,cbind)))
    tmp <- cbind(data.frame(1:dim(tmp)[1],object$historyfopt),tmp)
    names(tmp) <- c('Iteration','Cost Function',paste('x[',1:(dim(tmp)[2]-2),']',sep=''))
    cat('\nOptimization History:\n')
    print(tmp,row.names=FALSE)
  }
}
