print.DirichletRegData <- function(x, type=c("processed", "original"), ...){

  if(match.arg(type) == "processed") print(x[,]) else attr(x, "Y.original")

}
