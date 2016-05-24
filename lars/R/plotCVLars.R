plotCVLars <-
function(cv.lars.object,se=TRUE){
  mode=cv.lars.object$mode
  xlab=switch(mode,
    fraction="Fraction of final L1 norm",
    step="Number of steps"
    )
  index=cv.lars.object$index
  cv=cv.lars.object$cv
  cv.error=cv.lars.object$cv.error
      plot(index, cv, type = "b", ylim = range(cv, cv + cv.error, 
                                     cv - cv.error),xlab=xlab,ylab="Cross-Validated MSE")
    if(se)
      error.bars(index, cv + cv.error, cv - cv.error, 
                 width = 1/length(index))
invisible()
}


