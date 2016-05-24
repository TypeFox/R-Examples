print.transectHolder <- function(x, ...){
  print(paste("This object of class transectHolder contains ",
              length(x$transects), " transects.", sep=""))
  str(x)
}
