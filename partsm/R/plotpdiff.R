plotpdiff <- function(x){
  if (!(class(x) == "fit.piartsm"))
    stop("\n Object is not of class 'fit.piartsm'.\n")

  opar <- par(las=1)
  layout(matrix(c(1, 1, 2, 3), 2 , 2, byrow=TRUE))
  plot(x@pdiff.data, main="Periodically differenced data", ylab="", xlab="")
  bbplot(x@pdiff.data)
  monthplot(x@pdiff.data, ylab="")
  ##acf(x@pdiff.data, main="Autocorrelations", ylab="", na.action=na.pass)
  ##pacf(x@pdiff.data, main="Partial autocorrelations", ylab="", na.action=na.pass)
  par(opar)
}
