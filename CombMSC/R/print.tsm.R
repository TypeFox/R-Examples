`print.tsm` <-
function(x, ...){
  cat(paste(" "," ","AR", "I", "MA", sep='\t'))
  cat("\n")
  cat(paste( "Model:", " ",x$model[[1]], x$model[[2]], x$model[[3]], sep = '\t'))
  cat("\n")
  cat(paste("Seasonal:", x$seasonal$model[[1]], x$seasonal$model[[2]], x$seasonal$model[[3]], sep = '\t'))
  cat(paste("\n\nSeasonal Period:", x$seasonal$period, "\n", sep = '  '))
}

