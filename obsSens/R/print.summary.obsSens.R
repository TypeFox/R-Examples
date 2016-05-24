`print.summary.obsSens` <-
function(x,...){
  tmp <- x
  class(tmp) <- NULL
  attr(tmp,'log') <- NULL
  attr(tmp,'xname') <- NULL
  attr(tmp,'type') <- NULL
  if(attr(x,'log')){
    if( attr(x,'type') == 'cat' ){
      b <- 'on a Log Odds scale'
    } else if( attr(x, 'type') == 'num' ){
      b <- 'unscaled'
    } else if( attr(x, 'type') == 'surv' ){
      b <- 'on a Log Hazard Ratio scale'
    } else {
      warning('Unknown type')
      b <- ''
    }
  } else {
    if( attr(x,'type') == 'cat' ){
      b <- 'on an Odds Ratio scale'
    } else if( attr(x, 'type') == 'num' ){
      b <- 'On an Exponential scale'
    } else if( attr(x, 'type') == 'surv' ){
      b <- 'on a Hazard Ratio scale'
    } else {
      warning('Unknown type')
      b <- ''
    }
  }
  cat('Sensitivity analysis on variable', attr(x,'xname'),"\n",
      b, "\n\n")
  
  print.default(tmp, quote=FALSE,...)
  invisible(x)
}

