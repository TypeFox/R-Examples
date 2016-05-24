`plot.runJagsWrapper` <-
function(x, theoretical=FALSE, ...){
  
  if (class(x) != "runJagsWrapper")
    stop("'x' must be of class 'runJagsWrapper'")

  if (theoretical==TRUE) {
    a.plot <- plotFitted(x$seg.ratios, x$summary, theoretical=TRUE,
                         model=x$model, ...)
  } else {
    a.plot <- plotFitted(x$seg.ratios, x$summary, ...)
  }
  print(a.plot)
  return(a.plot)
}

